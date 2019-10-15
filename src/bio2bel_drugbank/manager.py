# -*- coding: utf-8 -*-

"""Defines the Bio2BEL DrugBank manager."""

import json
import logging
import os
import time
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Optional, TextIO, Tuple


import click
import networkx as nx
from sqlalchemy import func
from tqdm import tqdm

import bio2bel_hgnc
from bio2bel import AbstractManager
from bio2bel.manager.bel_manager import BELManagerMixin
from bio2bel.manager.flask_manager import FlaskMixin
from bio2bel.manager.namespace_manager import BELNamespaceManagerMixin
from pybel import BELGraph
from pybel.constants import ABUNDANCE, FUNCTION, IDENTIFIER, NAME, NAMESPACE, PROTEIN
from pybel.dsl import BaseEntity, abundance
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
    """Drug-target interactions."""

    module_name = MODULE_NAME
    _base = Base

    namespace_model = Drug
    identifiers_recommended = 'DrugBank'
    identifiers_pattern = r'^DB\d{5}$'
    identifiers_miriam = 'MIR:00000102'
    identifiers_namespace = 'drugbank'
    identifiers_url = 'http://identifiers.org/drugbank/'

    edge_model = DrugProteinInteraction

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

    def get_type_by_name(self, name: str) -> Optional[Type]:
        """Get a Type by name."""
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
        """List all chemical groups."""
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

    def get_or_create_species(self, name: str, tax_id: str) -> Species:
        """Get or create a Species by name."""
        m = self.species_to_model.get(name)
        if m is not None:
            return m

        m = self.get_species_by_name(name)
        if m is not None:
            self.species_to_model[name] = m
            return m

        m = self.species_to_model[name] = Species(name=name, tax_id=tax_id)
        self.session.add(m)
        return m

    def get_category_by_name(self, name: str) -> Optional[Category]:
        """Get a Category by name."""
        return self.session.query(Category).filter(Category.name == name).one_or_none()

    def get_or_create_category(self, name: str, **kwargs) -> Category:
        """Get or create a Category by name."""
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
        """Get or creates a Patent."""
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
        return 0 < self.count_drugs()

    def get_protein_by_uniprot_id(self, uniprot_id: str) -> Optional[Protein]:
        return self.session.query(Protein).filter(Protein.uniprot_id == uniprot_id).one_or_none()

    def get_protein_by_hgnc_id(self, hgnc_id: str) -> Optional[Protein]:
        res: List[Protein] = list(self.session.query(Protein).filter(Protein.hgnc_id == hgnc_id).all())

        if len(res) == 0:
            return
        if len(res) == 1:
            return res[0]

        warning_txt = '\n'.join(f' - {result.uniprot_id} {result.name}' for result in res)
        log.error('found multiple isoforms of hgnc:%s. Keeping first of:\n%s', hgnc_id, warning_txt)
        return res[0]

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

    def _create_drug_protein_interaction(self, drug: Drug, polypeptide, target):
        protein = self.get_or_create_protein(
            name=polypeptide['name'],
            uniprot_id=polypeptide['uniprot_id'],
            uniprot_accession=polypeptide['uniprot_accession'],
            species=self.get_or_create_species(
                name=polypeptide['organism'],
                tax_id=polypeptide['taxonomy'],
            ),
            hgnc_id=polypeptide.get('hgnc_id'),
            hgnc_symbol=polypeptide.get('hgnc_symbol'),
        )

        drug.protein_interactions.append(DrugProteinInteraction(
            drug=drug,
            protein=protein,
            known_action=(target['known_action'] == 'yes'),
            actions=[
                self.get_or_create_action(name.strip().lower())
                for name in target.get('actions', [])
            ],
            articles=[
                self.get_or_create_article(pubmed_id)
                for pubmed_id in target.get('articles', [])
            ],
            category=target['category'],
            type=target['type'],
        ))

    def populate(self, url: Optional[str] = None) -> None:
        """Populate DrugBank.

        :param url: Path to the DrugBank XML
        """
        root = get_xml_root(path=url)
        log.info(f'parsed DrugBank v{root.attrib["version"]}')

        # TODO get UniProt id to accession dictionary using Bio2BEL uniprot
        # TODO get HGNC identifier to symbol using Bio2BEL HGNC

        for drug_xml in tqdm(root, desc='Drugs'):
            drug_info = extract_drug_info(drug_xml)

            xrefs = drug_info['xrefs']
            pubchem_compound_id = None
            chebi_id = None
            for xref in xrefs:
                if xref['resource'] == 'PubChem Compound':
                    pubchem_compound_id = xref['identifier']
                elif xref['resource'] == 'ChEBI':
                    chebi_id = xref['identifier']

            drug = Drug(
                type=self.get_or_create_type(drug_info['type']),
                drugbank_id=drug_info['drugbank_id'],
                cas_number=drug_info['cas_number'],
                name=drug_info['name'],
                description=drug_info['description'],
                groups=[
                    self.get_or_create_group(name)
                    for name in drug_info['groups']
                ],
                atc_codes=[
                    AtcCode(name=name)
                    for name in drug_info['atc_codes']
                ],
                categories=[
                    self.get_or_create_category(**category)
                    for category in drug_info['categories']
                ],
                inchi=drug_info.get('inchi'),
                inchikey=drug_info.get('inchikey'),
                smiles=drug_info.get('smiles'),
                chebi_id=chebi_id,
                pubchem_compound_id=pubchem_compound_id,
                aliases=[
                    Alias(name=name)
                    for name in drug_info['aliases']
                ],
                patents=[
                    self.get_or_create_patent(**patent)
                    for patent in drug_info['patents']
                ],
                xrefs=[
                    DrugXref(**xref)
                    for xref in xrefs
                ],
            )

            drug.protein_interactions = []
            for row in drug_info['protein_interactions']:
                for polypeptide in row['polypeptides']:
                    self._create_drug_protein_interaction(drug, polypeptide, row['target'])

            self.session.add(drug)

        t = time.time()
        log.info('committing models')
        self.session.commit()
        log.info('committed models in %.2f seconds', time.time() - t)

    def count_drugs(self) -> int:
        """Count the number of drugs in the database."""
        return self._count_model(Drug)

    def list_drugs(self) -> List[Drug]:
        """List all drugs in the database."""
        return self._list_model(Drug)

    def count_types(self) -> int:
        """Count the number of types in the database."""
        return self._count_model(Type)

    def count_aliases(self) -> int:
        """Count the number of aliases in the database."""
        return self._count_model(Alias)

    def count_atc_codes(self) -> int:
        """Count the number of ATC codes in the database."""
        return self._count_model(AtcCode)

    def count_groups(self) -> int:
        """Count the number of groups in the database."""
        return self._count_model(Group)

    def count_categories(self) -> int:
        """Count the number of categories in the database."""
        return self._count_model(Category)

    def count_drugs_categories(self) -> int:
        """Count the number of drug-category relations in the database."""
        return self._count_model(drug_category)

    def count_drugs_groups(self) -> int:
        """Count the number of drug-group relations in the database."""
        return self._count_model(drug_group)

    def count_patents(self) -> int:
        """Count the number of patents in the database."""
        return self._count_model(Patent)

    def list_patents(self) -> List[Patent]:
        """List the patents in the database."""
        return self._list_model(Patent)

    def count_xrefs(self) -> int:
        """Count the number of cross-references in the database."""
        return self._count_model(DrugXref)

    def get_xrefs_by_resource(self, resource) -> List[DrugXref]:
        return self.session.query(DrugXref).filter(DrugXref.resource == resource).all()

    def summarize_xrefs(self) -> Counter:
        return Counter(
            self.session.query(DrugXref.resource, func.count(DrugXref.resource)).group_by(DrugXref.resource).all())

    def count_species(self) -> int:
        return self._count_model(Species)

    def count_proteins(self) -> int:
        return self._count_model(Protein)

    def list_proteins(self) -> List[Protein]:
        return self._list_model(Protein)

    def count_actions(self) -> int:
        return self._count_model(Action)

    def count_drug_protein_interactions(self) -> int:
        return self._count_model(DrugProteinInteraction)

    def list_drug_protein_interactions(self) -> List[DrugProteinInteraction]:
        """List drug-protein interactions."""
        return self._list_model(DrugProteinInteraction)

    def count_articles(self) -> int:
        return self._count_model(Article)

    def summarize(self) -> Dict[str, int]:
        """Summarize the database."""
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
    def _get_identifier(drug: Drug) -> str:
        return drug.drugbank_id

    def _create_namespace_entry_from_model(self, model: Drug, namespace: Namespace) -> NamespaceEntry:
        return NamespaceEntry(
            encoding='A',
            name=model.name,
            identifier=model.drugbank_id,
            namespace=namespace,
        )

    def lookup_target(self, node: BaseEntity) -> Optional[Protein]:
        namespace = node.get(NAMESPACE)
        if node[FUNCTION] != PROTEIN or namespace is None:
            return

        identifier = node.get(IDENTIFIER)
        if namespace.lower() == 'hgnc' and identifier:
            return self.get_protein_by_hgnc_id(identifier)

        if namespace.lower() == 'uniprot' and identifier:
            return self.get_protein_by_uniprot_id(identifier)

    def iter_targets(self, graph: BELGraph) -> Iterable[Tuple[BaseEntity, Protein]]:
        for node in graph:
            protein_model = self.lookup_target(node)
            if protein_model is not None:
                yield node, protein_model

    def enrich_targets(self, graph: BELGraph) -> None:
        """Enrich the protein targets in the graph with Drug-Protein interactions from DrugBank."""
        self.add_namespace_to_graph(graph)

        c = 0
        for node_data, protein_model in list(self.iter_targets(graph)):
            for dpi in protein_model.drug_interactions:
                dpi._add_to_graph(graph, dpi.drug.as_bel(), node_data)
                c += 1

        log.info('added %d drug-protein interactions.', c)

    def get_drug_by_drugbank_id(self, drugbank_id: str) -> Optional[Drug]:
        return self.session.query(Drug).filter(Drug.drugbank_id == drugbank_id).one_or_none()

    def get_drug_by_name(self, name) -> Optional[Drug]:
        return self.session.query(Drug).filter(Drug.name == name).one_or_none()

    def get_drug_by_inchi(self, inchi: str) -> Optional[Drug]:
        return self.session.query(Drug).filter(Drug.inchi == inchi).one_or_none()

    def get_drug_by_inchikey(self, inchikey: str) -> Optional[Drug]:
        return self.session.query(Drug).filter(Drug.inchikey == inchikey).one_or_none()

    def get_drug_by_xref(self, resource: str, identifier: str) -> Optional[Drug]:
        # need to join the xref table for this one
        xref = self.session.query(DrugXref).filter(DrugXref.has_identifier(resource, identifier)).one_or_none()
        if xref:
            return xref.drug

    def lookup_drug(self, node: BaseEntity) -> Optional[Drug]:
        """Try and look up a drug."""
        namespace = node.get(NAMESPACE)

        if node[FUNCTION] != ABUNDANCE or namespace is None:
            return

        name, identifier = node.get(NAME), node.get(IDENTIFIER)

        if namespace.lower() == 'drugbank':
            if identifier is not None:
                return self.get_drug_by_drugbank_id(identifier)
            elif name.startswith('DB'):
                return self.get_drug_by_drugbank_id(name)
            else:
                return self.get_drug_by_name(name)

    def iter_drugs(self, graph: BELGraph, use_tqdm: bool = False) -> Iterable[Tuple[BaseEntity, Drug]]:
        """Iterate over the drugs in the graph."""
        it = (
            tqdm(graph, desc='DrugBank chemicals')
            if use_tqdm else
            graph
        )
        for node in it:
            drug_model = self.lookup_drug(node)
            if drug_model is not None:
                yield node, drug_model

    def normalize_drugs(self, graph: BELGraph, use_tqdm: bool = False) -> None:
        """Normalize the drugs in the graph."""
        mapping = {
            node: drug_model.as_bel()
            for node, drug_model in self.iter_drugs(graph, use_tqdm=use_tqdm)
        }
        nx.relabel_nodes(graph, mapping, copy=False)

    def enrich_drug_inchi(self, graph: BELGraph) -> None:
        """Enrich drugs in the graph with their InChI equivalent nodes."""
        self.add_namespace_to_graph(graph)

        for node, drug_model in list(self.iter_drugs(graph)):
            if drug_model.inchi:
                graph.add_equivalence(node, drug_model.as_inchi_bel())

    def enrich_drug_equivalences(self, graph: BELGraph) -> None:
        """Enrich drugs in the graph with their equivalent nodes."""
        self.add_namespace_to_graph(graph)

        for node, drug_model in list(self.iter_drugs(graph)):
            if drug_model.inchi:
                graph.add_equivalence(node, drug_model.as_inchi_bel())

            if drug_model.inchikey:
                graph.add_equivalence(node, drug_model.as_inchikey_bel())

            for xref in drug_model.xrefs:
                resource = xref.resource.lower()
                identifier = xref.identifier

                if xref.resource in {'chebi', 'chembl'}:
                    graph.add_equivalence(node, abundance(namespace=resource, identifier=identifier))
                elif xref.resource == 'KEGG Compound':
                    # https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000013
                    graph.add_equivalence(node, abundance(namespace='kegg.compound', identifier=identifier))
                elif xref.resource == 'PubChem Substance':
                    # https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000033
                    graph.add_equivalence(node, abundance(namespace='pubchem.substance', identifier=identifier))
                elif xref.resource == 'PubChem Compound':
                    # https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000034
                    graph.add_equivalence(node, abundance(namespace='pubchem.compound', identifier=identifier))

                # TODO there are plenty more. implement as other bio2bel repositories need

    def enrich_drugs(self, graph: BELGraph) -> None:
        """Enrich drugs in the graph with their targets."""
        self.add_namespace_to_graph(graph)

        for node_data, drug_model in list(self.iter_drugs(graph)):
            for dpi in drug_model.protein_interactions:
                dpi.add_to_graph(graph)

    def to_bel(self, drug_namespace: Optional[str] = None, target_namespace: Optional[str] = None) -> BELGraph:
        """Export DrugBank as BEL."""
        graph = BELGraph(
            name='DrugBank',
            version='5.1.4',
        )

        self.add_namespace_to_graph(graph)

        hgnc_manager = bio2bel_hgnc.Manager(engine=self.engine, session=self.session)
        hgnc_manager.add_namespace_to_graph(graph)

        dpis = self.list_drug_protein_interactions()
        dpis: Iterable[DrugProteinInteraction] = tqdm(
            dpis,
            total=self.count_drug_protein_interactions(),
            desc='Mapping drug-protein interactions to BEL',
        )
        for dpi in dpis:
            dpi.add_to_graph(graph, drug_namespace=drug_namespace, target_namespace=target_namespace)

        return graph

    def get_hgnc_id_to_drugs(self) -> Dict[str, List[str]]:
        """Get a dictionary of HGNC identifiers (not prepended with HGNC:) to list of drug names."""
        rv = defaultdict(list)

        for dpi in tqdm(self.list_drug_protein_interactions(),
                        total=self.count_drug_protein_interactions(),
                        desc='getting DTIs'):

            if dpi.protein.hgnc_id is None:
                continue

            drug_name = dpi.drug.name
            hgnc_id = dpi.protein.hgnc_id[len('HGNC:'):]

            rv[hgnc_id].append(drug_name)

        return dict(rv)

    def get_drug_to_hgnc_ids(self, cache=True, recalculate=False) -> Dict[str, List[str]]:
        """Get a dictionary of drug names to lists HGNC identifiers (not prepended with HGNC:)."""
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

        return dict(rv)

    def get_drug_to_hgnc_symbols(self, cache=True, recalculate=False) -> Dict[str, List[str]]:
        """Get a dictionary of drug names to HGNC gene symbols."""
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

        return dict(rv)

    def get_interactions_by_hgnc_id(self, hgnc_id: str) -> List[DrugProteinInteraction]:
        """Get the drug targets for a given HGNC identifier."""
        protein = self.get_protein_by_hgnc_id(hgnc_id)

        if not protein:
            return []

        return [
            interaction
            for interaction in protein.drug_interactions
        ]

    @classmethod
    def get_cli(cls) -> click.Group:
        """Append the lister."""
        main = super().get_cli()

        @main.command()
        @click.option('-f', '--file', type=click.File('w'))
        @click.pass_obj
        def export(manager: Manager, file: TextIO):
            """Export DTIs as a tall/skinny."""
            print(
                'drug_drugbank_id',
                'drug_drugbank_name',
                'drug_type',
                'drug_drugbank_inchikey',
                'target_entrez_id',
                'target_hgnc_id',
                'target_uniprot_id',
                'target_uniprot_accession',
                'target_name',
                'category',
                sep='\t',
                file=file,
            )
            for dpi in tqdm(manager.list_drug_protein_interactions(), total=manager.count_drug_protein_interactions()):
                print(
                    dpi.drug.drugbank_id,
                    dpi.drug.name,
                    dpi.drug.type.name,
                    dpi.drug.inchikey,
                    dpi.protein.entrez_id,
                    dpi.protein.hgnc_id,
                    dpi.protein.uniprot_id,
                    dpi.protein.uniprot_accession,
                    dpi.protein.name,
                    dpi.category,
                    sep='\t',
                    file=file,
                )

        return main
