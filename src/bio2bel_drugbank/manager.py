# -*- coding: utf-8 -*-

import logging

import time
from tqdm import tqdm

from bio2bel import AbstractManager
from pybel.manager.models import Namespace, NamespaceEntry
from .constants import MODULE_NAME
from .models import Alias, AtcCode, Base, Category, Drug, Group, Patent, Type, Xref, drug_category, drug_group
from .parser import extract_drug_info, get_xml_root

__all__ = ['Manager']

log = logging.getLogger(__name__)


class Manager(AbstractManager):
    """Manager for Bio2BEL DrugBank"""
    module_name = MODULE_NAME
    flask_admin_models = [Drug, Alias, AtcCode, Category, Group, Type, Patent, Xref]

    def __init__(self, connection=None):
        super().__init__(connection=connection)

        self.type_to_model = {}
        self.group_to_model = {}
        self.category_to_model = {}
        self.patent_to_model = {}

    @property
    def base(self):
        return Base

    def get_type_by_name(self, name):
        return self.session.query(Type).filter(Type.name == name).one_or_none()

    def get_or_create_type(self, name):
        m = self.type_to_model.get(name)
        if m is not None:
            return m

        m = self.get_type_by_name(name)
        if m is not None:
            self.type_to_model[name] = m
            return m

        m = self.type_to_model[name] = Type(name=name)
        self.session.add(m)
        return m

    def get_group_by_name(self, name):
        return self.session.query(Group).filter(Group.name == name).one_or_none()

    def get_or_create_group(self, name):
        m = self.group_to_model.get(name)
        if m is not None:
            return m

        m = self.get_group_by_name(name)
        if m is not None:
            self.group_to_model[name] = m
            return m

        m = self.group_to_model[name] = Group(name=name)
        self.session.add(m)
        return m

    def get_category_by_name(self, name):
        return self.session.query(Category).filter(Category.name == name).one_or_none()

    def get_or_create_category(self, name, **kwargs):
        m = self.category_to_model.get(name)
        if m is not None:
            return m

        m = self.get_category_by_name(name)
        if m is not None:
            self.category_to_model[name] = m
            return m

        m = self.category_to_model[name] = Category(name=name, **kwargs)
        self.session.add(m)
        return m

    def get_or_create_patent(self, country, patent_id, **kwargs):
        m = self.patent_to_model.get((country, patent_id))
        if m is not None:
            return m

        m = self.session.query(Patent).filter(Patent.filter_pk(country, patent_id)).one_or_none()
        if m is not None:
            self.patent_to_model[(country, patent_id)] = m
            return m

        m = self.patent_to_model[(country, patent_id)] = Patent(
            country=country,
            patent_id=patent_id,
            **kwargs
        )
        self.session.add(m)
        return m

    def is_populated(self):
        """Checks if the databse is populated by counting the drugs

        :rtype: bool
        """
        return 0 != self.count_drugs()

    def populate(self, url=None):
        """Populates DrugBank

        :param Optional[str] url: Path to the DrugBank XML
        """
        root = get_xml_root(url=url)

        log.info('building models')

        for drug_xml in tqdm(root):
            drug = extract_drug_info(drug_xml)

            drug_model = Drug(
                type=self.get_or_create_type(drug['type']),
                drugbank_id=drug['drugbank_id'],
                cas_number=drug['cas_number'],
                name=drug['name'],
                description=drug['description'],
                groups=[
                    self.get_or_create_group(name)
                    for name in drug['groups']
                ],
                atc_codes=[
                    AtcCode(name=name)
                    for name in drug['atc_codes']
                ],
                categories=[
                    self.get_or_create_category(**category)
                    for category in drug['categories']
                ],
                inchi=drug.get('inchi'),
                inchikey=drug.get('inchikey'),
                aliases=[
                    Alias(name=name)
                    for name in drug['aliases']
                ],
                patents=[
                    self.get_or_create_patent(**patent)
                    for patent in drug['patents']
                ],
                xrefs=[
                    Xref(**xref)
                    for xref in drug['xrefs']
                ]
            )

            self.session.add(drug_model)

        t = time.time()
        log.info('committing models')
        self.session.commit()
        log.info('committed models in %.2f seconds', time.time() - t)

    def count_drugs(self):
        """Count the number of drugs in the database

        :rtype: int
        """
        return self._count_model(Drug)

    def count_types(self):
        """Count the number of types in the database

        :rtype: int
        """
        return self._count_model(Type)

    def count_alises(self):
        """Count the number of aliases in the database

        :rtype: int
        """
        return self._count_model(Alias)

    def count_atc_codes(self):
        """Count the number of ATC codes in the database

        :rtype: int
        """
        return self._count_model(AtcCode)

    def count_groups(self):
        """Count the number of groups in the database

        :rtype: int
        """
        return self._count_model(Group)

    def count_categories(self):
        """Count the number of categories in the database

        :rtype: int
        """
        return self._count_model(Category)

    def count_drugs_categories(self):
        """Count the number of drug-category relations in the database

        :rtype: int
        """
        return self._count_model(drug_category)

    def count_drugs_groups(self):
        """Count the number of drug-group relations in the database

        :rtype: int
        """
        return self._count_model(drug_group)

    def count_patents(self):
        """Count the number of patents in the database

        :rtype: int
        """
        return self._count_model(Patent)

    def list_patents(self):
        """Lists the patents in the database

        :rtype: list[Patent]
        """
        return self._list_model(Patent)

    def count_xrefs(self):
        """Count the number of cross-references in the database

        :rtype: int
        """
        return self._count_model(Xref)

    def summarize(self):
        """Summarizes the database

        :rtype: dict[str,int]
        """
        return dict(
            drugs=self.count_drugs(),
            types=self.count_types(),
            aliases=self.count_alises(),
            atc_codes=self.count_atc_codes(),
            groups=self.count_groups(),
            categories=self.count_categories(),
            patents=self.count_patents(),
            xrefs=self.count_xrefs(),
        )

    def _iterate_id_name(self):
        return tqdm(self.session.query(Drug.drugbank_id, Drug.name), total=self.count_drugs())

    def _get_namespace_entries(self):
        return [
            NamespaceEntry(encoding='A', identifier=identifier, name=name)
            for identifier, name in self._iterate_id_name()
        ]

    def _make_namespace(self):
        """
        :rtype: pybel.manager.models.Namespace
        """
        entries = self._get_namespace_entries()
        _namespace_keyword = self._get_namespace_keyword()
        ns = Namespace(
            name=_namespace_keyword,
            keyword=_namespace_keyword,
            url=_namespace_keyword,
            version=str(time.asctime()),
            entries=entries,
        )
        self.session.add(ns)

        t = time.time()
        log.info('committing models')
        self.session.commit()
        log.info('committed models in %.2f seconds', time.time() - t)

        return ns

    def _update_namespace(self, ns):
        """
        :param pybel.manager.models.Namespace ns:
        """
        old = {term.identifier for term in ns.entries}
        new_count = 0

        for identifier, name in self._iterate_id_name():
            if identifier in old:
                continue

            new_count += 1
            entry = NamespaceEntry(encoding='A', identifier=identifier, name=name, namespace=ns)
            self.session.add(entry)

        t = time.time()
        log.info('got %d new entries. committing models', new_count)
        self.session.commit()
        log.info('committed models in %.2f seconds', time.time() - t)

    def upload_bel_namespace(self):
        """
        :rtype: pybel.manager.models.Namespace
        """
        if not self.is_populated():
            self.populate()

        ns = self._get_default_namespace()

        if ns is None:
            return self._make_namespace()

        self._update_namespace(ns)

        return ns
