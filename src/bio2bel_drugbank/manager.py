# -*- coding: utf-8 -*-

"""Defines the Bio2BEL DrugBank manager."""

import json
import logging
import os
from collections import defaultdict

import bio2bel_hgnc
import time
from sqlalchemy import func
from tqdm import tqdm
from typing import List, Optional

from bio2bel import AbstractManager
from bio2bel.manager.bel_manager import BELManagerMixin
from bio2bel.manager.flask_manager import FlaskMixin
from bio2bel.manager.namespace_manager import BELNamespaceManagerMixin
from pybel import BELGraph
from pybel.manager.models import Namespace, NamespaceEntry
from .constants import DATA_DIR, MODULE_NAME
from .models import (
    Action, Alias, Article, AtcCode, Base, Category, Drug, DrugProteinInteraction, DrugXref, Group, Patent, Protein,
    Species, Type, drug_category, drug_group,
)
from .parser import extract_drug_info, get_xml_root

__all__ = ['Manager']

log = logging.getLogger(__name__)

_dti_ids_cache_path = os.path.join(DATA_DIR, 'drug_to_gene_ids.json')
_dti_symbols_cache_path = os.path.join(DATA_DIR, 'drug_to_gene_symbols.json')


class Manager(AbstractManager, FlaskMixin, BELManagerMixin, BELNamespaceManagerMixin):
    """Manager for Bio2BEL DrugBank."""

    module_name = MODULE_NAME

    namespace_model = Drug
    identifiers_recommended = 'DrugBank'
    identifiers_pattern = '^DB\d{5}$'
    identifiers_miriam = 'MIR:00000102'
    identifiers_namespace = 'drugbank'
    identifiers_url = 'http://identifiers.org/drugbank/'

    flask_admin_models = [Drug, Alias, AtcCode, Category, Group, Type, Patent, DrugXref, Species, Protein,
                          DrugProteinInteraction, Action, Article]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.type_to_model = {}
        self.group_to_model = {}
        self.category_to_model = {}
        self.patent_to_model = {}
        self.species_to_model = {}
        self.action_to_model = {}
        self.pmid_to_model = {}
        self.uniprot_id_to_protein = {}

    @property
    def _base(self):
        return Base

    def get_type_by_name(self, name: str) -> Optional[Type]:
        """Get a Type by name.

        :param str name:
        """
        return self.session.query(Type).filter(Type.name == name).one_or_none()

    def get_or_create_type(self, name: str) -> Type:
        """Get or create a Type by name."""
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

    def list_groups(self) -> List[Group]:
        return self._list_model(Group)

    def get_group_by_name(self, name: str) -> Optional[Group]:
        """Get a Group by name."""
        return self.session.query(Group).filter(Group.name == name).one_or_none()

    def get_or_create_group(self, name: str) -> Group:
        """Get or create a Group by name."""
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

    def get_species_by_name(self, name: str) -> Optional[Species]:
        """Get a Species by name."""
        return self.session.query(Species).filter(Species.name == name).one_or_none()

    def get_or_create_species(self, name: str) -> Species:
        """Get or create a Species by name."""
        m = self.species_to_model.get(name)
        if m is not None:
            return m

        m = self.get_species_by_name(name)
        if m is not None:
            self.species_to_model[name] = m
            return m

        m = self.species_to_model[name] = Species(name=name)
        self.session.add(m)
        return m

    def get_category_by_name(self, name: str) -> Optional[Category]:
        """Get a Category by name."""
        return self.session.query(Category).filter(Category.name == name).one_or_none()

    def get_or_create_category(self, name: str, **kwargs) -> Category:
        """Gets or creates a Category by name."""
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

    def get_or_create_patent(self, country: str, patent_id: str, **kwargs) -> Patent:
        """Gets or creates a Patent."""
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

    def is_populated(self) -> bool:
        """Check if the database is populated by counting the drugs."""
        return 0 != self.count_drugs()

    def get_protein_by_uniprot_id(self, uniprot_id: str) -> Optional[Protein]:
        return self.session.query(Protein).filter(Protein.uniprot_id == uniprot_id).one_or_none()

    def get_protein_by_hgnc_id(self, hgnc_id):
        return self.session.query(Protein).filter(Protein.hgnc_id == hgnc_id).one_or_none()

    def get_or_create_protein(self, uniprot_id: str, **kwargs) -> Protein:
        m = self.uniprot_id_to_protein.get(uniprot_id)
        if m is not None:
            return m

        m = self.get_protein_by_uniprot_id(uniprot_id)
        if m is not None:
            self.uniprot_id_to_protein[uniprot_id] = m
            return m

        m = self.uniprot_id_to_protein[uniprot_id] = Protein(
            uniprot_id=uniprot_id,
            **kwargs
        )
        self.session.add(m)
        return m

    def get_action_by_name(self, name: str) -> Optional[Action]:
        return self.session.query(Action).filter(Action.name == name).one_or_none()

    def get_or_create_action(self, name: str) -> Action:
        m = self.action_to_model.get(name)
        if m is not None:
            return m

        m = self.get_action_by_name(name)
        if m is not None:
            self.action_to_model[name] = m
            return m

        m = self.action_to_model[name] = Action(name=name)
        self.session.add(m)
        return m

    def get_article_by_pmid(self, pubmed_id: str):
        return self.session.query(Article).filter(Article.pubmed_id == pubmed_id).one_or_none()

    def get_or_create_article(self, pubmed_id) -> Article:
        m = self.pmid_to_model.get(pubmed_id)
        if m is not None:
            return m

        m = self.get_article_by_pmid(pubmed_id)
        if m is not None:
            self.pmid_to_model[pubmed_id] = m
            return m

        m = self.pmid_to_model[pubmed_id] = Article(pubmed_id=pubmed_id)
        self.session.add(m)
        return m

    def _create_drug_protein_interaction(self, drug_model, d):
        """

        :param dict d:
        :return: DrugProteinInteraction
        """
        protein = self.get_or_create_protein(
            uniprot_id=d['uniprot_id'],
            species=self.get_or_create_species(d['organism']),
            name=d.get('name'),
            hgnc_id=d.get('hgnc_id')
        )

        dpi = DrugProteinInteraction(
            drug=drug_model,
            protein=protein,
            known_action=(d['known_action'] == 'yes'),
            actions=[self.get_or_create_action(name.strip().lower()) for name in d.get('actions', [])],
            articles=[self.get_or_create_article(pubmed_id) for pubmed_id in d.get('articles', [])],
            category=d['category'],
        )
        self.session.add(dpi)
        return dpi

    def populate(self, url=None):
        """Populates DrugBank

        :param Optional[str] url: Path to the DrugBank XML
        """
        root = get_xml_root(url=url)

        log.info('building models')

        for drug_xml in tqdm(root, desc='Drugs'):
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
                    DrugXref(**xref)
                    for xref in drug['xrefs']
                ]
            )

            drug_model.protein_interactions = [
                self._create_drug_protein_interaction(drug_model, drug_protein_interaction_data)
                for drug_protein_interaction_data in drug['protein_interactions']
            ]

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

    def list_drugs(self):
        """List all drugs in the database

        :rtype: list[Drug]
        """
        return self._list_model(Drug)

    def count_types(self):
        """Count the number of types in the database

        :rtype: int
        """
        return self._count_model(Type)

    def count_aliases(self):
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

    def count_drugs_categories(self) -> int:
        """Count the number of drug-category relations in the database

        :rtype: int
        """
        return self._count_model(drug_category)

    def count_drugs_groups(self) -> int:
        """Count the number of drug-group relations in the database

        :rtype: int
        """
        return self._count_model(drug_group)

    def count_patents(self) -> int:
        """Count the number of patents in the database

        :rtype: int
        """
        return self._count_model(Patent)

    def list_patents(self):
        """Lists the patents in the database

        :rtype: list[Patent]
        """
        return self._list_model(Patent)

    def count_xrefs(self) -> int:
        """Count the number of cross-references in the database

        :rtype: int
        """
        return self._count_model(DrugXref)

    def get_xrefs_by_resource(self, resource) -> List[DrugXref]:
        return self.session.query(DrugXref).filter(DrugXref.resource == resource).all()

    def summarize_xrefs(self):
        return self.session.query(DrugXref.resource, func.count(DrugXref.resource)).group_by(DrugXref.resource).all()

    def count_species(self) -> int:
        return self._count_model(Species)

    def count_proteins(self) -> int:
        return self._count_model(Protein)

    def list_proteins(self) -> List[Protein]:
        return self._list_model(Protein)

    def count_actions(self) -> int:
        return self._count_model(Action)

    def count_drug_protein_interactions(self):
        return self._count_model(DrugProteinInteraction)

    def list_drug_protein_interactions(self) -> List[DrugProteinInteraction]:
        """List drug-protein interactions

        :rtype: list[DrugProteinInteraction]
        """
        return self._list_model(DrugProteinInteraction)

    def count_articles(self) -> int:
        return self._count_model(Article)

    def summarize(self):
        """Summarizes the database

        :rtype: dict[str,int]
        """
        return dict(
            drugs=self.count_drugs(),
            types=self.count_types(),
            aliases=self.count_aliases(),
            atc_codes=self.count_atc_codes(),
            groups=self.count_groups(),
            categories=self.count_categories(),
            patents=self.count_patents(),
            xrefs=self.count_xrefs(),
            proteins=self.count_proteins(),
            species=self.count_species(),
            actions=self.count_actions(),
            drug_protein_interactions=self.count_drug_protein_interactions(),
        )

    @staticmethod
    def _get_identifier(drug: Drug):
        return drug.drugbank_id

    def _create_namespace_entry_from_model(self, model: Drug, namespace: Namespace):
        return NamespaceEntry(encoding='A', name=model.name, identifier=model.drugbank_id, namespace=namespace)

    def to_bel(self):
        graph = BELGraph(
            name='DrugBank',
            version='5.1',
        )

        drugbank_namespace = self.upload_bel_namespace()
        graph.namespace_url[drugbank_namespace.keyword] = drugbank_namespace.url

        hgnc_manager = bio2bel_hgnc.Manager(engine=self.engine, session=self.session)
        hgnc_namespace = hgnc_manager.upload_bel_namespace()
        graph.namespace_url[hgnc_namespace.keyword] = hgnc_namespace.url

        dpis = self.list_drug_protein_interactions()
        for dpi in tqdm(dpis, total=self.count_drug_protein_interactions(),
                        desc='Mapping drug-protein interactions to BEL'):
            if dpi.protein.hgnc_id is None:
                continue
            dpi.add_to_graph(graph)

        return graph

    def get_hgnc_id_to_drugs(self):
        """Gets a dictionary of HGNC identifiers (not prepended with HGNC:) to list of drug names

        :rtype: dict[str,list[str]]
        """
        rv = defaultdict(list)

        for dpi in tqdm(self.list_drug_protein_interactions(),
                        total=self.count_drug_protein_interactions(),
                        desc='getting DTIs'):

            if dpi.protein.hgnc_id is None:
                continue

            drug_name = dpi.drug.name
            hgnc_id = dpi.protein.hgnc_id[len('HGNC:'):]

            rv[hgnc_id].append(drug_name)

        return rv

    def get_drug_to_hgnc_ids(self, cache=True, recalculate=False):
        """Gets a dictionary of drug names to lists HGNC identifiers (not prepended with HGNC:)

        :rtype: dict[str,list[str]]
        """
        if cache and not recalculate and os.path.exists(_dti_ids_cache_path):
            log.info('loading cached DTIs')
            with open(_dti_ids_cache_path) as file:
                return json.load(file)

        rv = defaultdict(list)

        for dpi in tqdm(self.list_drug_protein_interactions(),
                        total=self.count_drug_protein_interactions(),
                        desc='getting DTIs'):

            if dpi.protein.hgnc_id is None:
                continue

            drug_name = dpi.drug.name
            hgnc_id = dpi.protein.hgnc_id[len('HGNC:'):]

            rv[drug_name].append(hgnc_id)

        if cache:
            with open(_dti_ids_cache_path, 'w') as file:
                log.info('dumping cached DTIs')
                json.dump(rv, file)

        return rv

    def get_drug_to_hgnc_symbols(self, cache=True, recalculate=False):
        """Gets a dictionary of drug names to HGNC gene symbols.

        Requires the installation of ``bio2bel_hgnc``

        :rtype: dict[str,list[str]]
        :raises: ImportError
        """
        if cache and not recalculate and os.path.exists(_dti_symbols_cache_path):
            log.debug('loading cached DTIs with gene symbols')
            with open(_dti_symbols_cache_path) as file:
                return json.load(file)

        hgnc_manager = bio2bel_hgnc.Manager(engine=self.engine, session=self.session)
        if not hgnc_manager.is_populated():
            hgnc_manager.populate()

        hgnc_id_symbol_mapping = hgnc_manager.build_hgnc_id_symbol_mapping()
        drug_to_hgnc_ids = self.get_drug_to_hgnc_ids()

        rv = defaultdict(list)

        for drug, hgnc_ids in drug_to_hgnc_ids.items():
            for hgnc_id in hgnc_ids:
                hgnc_symbol = hgnc_id_symbol_mapping.get(hgnc_id)

                if hgnc_symbol is None:
                    log.warning('could not map HGNC identifier: %s', hgnc_id)
                    continue

                rv[drug].append(hgnc_symbol)

        if cache:
            with open(_dti_symbols_cache_path, 'w') as file:
                log.info('dumping cached DTIs')
                json.dump(rv, file)

        return rv

    def get_interactions_by_hgnc_id(self, hgnc_id):
        """Gets the drug targets for a given HGNC identifier.

        :param str hgnc_id: HGNC identifier
        :rtype: Optional[list[DrugProteinInteraction]]
        """
        protein = self.get_protein_by_hgnc_id(hgnc_id)

        if not protein:
            return None

        return [
            interaction
            for interaction in protein.drug_interactions
        ]

    # TODO: ADD enrich graph methods
