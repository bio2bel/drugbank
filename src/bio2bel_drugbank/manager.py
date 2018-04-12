# -*- coding: utf-8 -*-

from bio2bel import AbstractManager
from .constants import MODULE_NAME
from .models import Alias, AtcCode, Base, Category, Drug, Group, Patent, Type, drug_category, drug_group, Xref
from .parser import *
import time

__all__ = ['Manager']

log = logging.getLogger(__name__)


class Manager(AbstractManager):
    module_name = MODULE_NAME
    flask_admin_models = [Drug, Alias, AtcCode, Category, Group, Type, Patent]

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

    def populate(self, url=None):
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

    def count_xrefs(self):
        """Count the number of cross-references in the database

        :rtype: int
        """
        return self._count_model(Xref)

    def summarize(self):
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
