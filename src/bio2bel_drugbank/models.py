# -*- coding: utf-8 -*-

"""Database models for bio2bel_drugbank."""

from typing import Set

from sqlalchemy import Boolean, Column, Date, ForeignKey, Integer, String, Table, Text, and_
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import backref, relationship

from pybel import BELGraph
from pybel.constants import REGULATES
from pybel.dsl import BaseEntity, abundance, activity, protein
from .constants import MODULE_NAME
from .patent_utils import download_google_patents

__all__ = [
    'Action',
    'Base',
    'Type',
    'Drug',
    'DrugProteinInteraction',
    'Alias',
    'AtcCode',
    'Group',
    'Category',
    'Patent',
    'DrugXref',
    'Protein',
    'Species',
    'Article',
    'drug_category',
    'drug_group',
    'drug_patent',
]

Base = declarative_base()

DRUG_TABLE_NAME = f'{MODULE_NAME}_drug'
TYPE_TABLE_NAME = f'{MODULE_NAME}_type'
ALIAS_TABLE_NAME = f'{MODULE_NAME}_alias'
ATC_TABLE_NAME = f'{MODULE_NAME}_atc'
CATEGORY_TABLE_NAME = f'{MODULE_NAME}_category'
GROUP_TABLE_NAME = f'{MODULE_NAME}_group'
DRUG_CATEGORY_TABLE_NAME = f'{MODULE_NAME}_drug_category'
DRUG_GROUP_TABLE_NAME = f'{MODULE_NAME}_drug_group'
PATENT_TABLE_NAME = f'{MODULE_NAME}_patent'
DRUG_PATENT_TABLE_NAME = f'{MODULE_NAME}_drug_patent'
DRUG_XREF_TABLE_NAME = f'{MODULE_NAME}_drug_xref'
SPECIES_TABLE_NAME = f'{MODULE_NAME}_species'
PROTEIN_TABLE_NAME = f'{MODULE_NAME}_protein'
ACTION_TABLE_NAME = f'{MODULE_NAME}_action'
DRUG_PROTEIN_TABLE_NAME = f'{MODULE_NAME}_drug_protein'
DRUG_PROTEIN_ACTION_TABLE_NAME = f'{MODULE_NAME}_drug_protein_action'
ARTICLE_TABLE_NAME = f'{MODULE_NAME}_article'
DRUG_PROTEIN_ARTICLE_TABLE_NAME = f'{MODULE_NAME}_drug_protein_article'

drug_group = Table(
    DRUG_GROUP_TABLE_NAME,
    Base.metadata,
    Column('drug_id', Integer, ForeignKey(f'{DRUG_TABLE_NAME}.id'), primary_key=True),
    Column('group_id', Integer, ForeignKey(f'{GROUP_TABLE_NAME}.id'), primary_key=True)
)

drug_category = Table(
    DRUG_CATEGORY_TABLE_NAME,
    Base.metadata,
    Column('drug_id', Integer, ForeignKey(f'{DRUG_TABLE_NAME}.id'), primary_key=True),
    Column('category_id', Integer, ForeignKey(f'{CATEGORY_TABLE_NAME}.id'), primary_key=True)
)

drug_patent = Table(
    DRUG_PATENT_TABLE_NAME,
    Base.metadata,
    Column('drug_id', Integer, ForeignKey(f'{DRUG_TABLE_NAME}.id'), primary_key=True),
    Column('patent_id', Integer, ForeignKey(f'{PATENT_TABLE_NAME}.id'), primary_key=True)
)

dpi_action = Table(
    DRUG_PROTEIN_ACTION_TABLE_NAME,
    Base.metadata,
    Column('drug_protein_id', Integer, ForeignKey(f'{DRUG_PROTEIN_TABLE_NAME}.id'), primary_key=True),
    Column('action_id', Integer, ForeignKey(f'{ACTION_TABLE_NAME}.id'), primary_key=True)
)

dpi_article = Table(
    DRUG_PROTEIN_ARTICLE_TABLE_NAME,
    Base.metadata,
    Column('drug_protein_id', Integer, ForeignKey(f'{DRUG_PROTEIN_TABLE_NAME}.id'), primary_key=True),
    Column('article_id', Integer, ForeignKey(f'{ARTICLE_TABLE_NAME}.id'), primary_key=True)
)

PATENT_PREFIX_MAP = {
    'United States': 'US',
    'Canada': 'CA',
}


class Type(Base):
    """Represents the type of a drug - either small molecule or biologic."""

    __tablename__ = TYPE_TABLE_NAME
    id = Column(Integer, primary_key=True)

    name = Column(String(255), unique=True, index=True, nullable=False)

    def __repr__(self):
        return self.name


class Drug(Base):
    """Represents a drug."""

    __tablename__ = DRUG_TABLE_NAME
    id = Column(Integer, primary_key=True)

    type_id = Column(Integer, ForeignKey(f'{TYPE_TABLE_NAME}.id'), nullable=False)
    type = relationship(Type)

    drugbank_id = Column(String(255), nullable=False)
    cas_number = Column(String(255), nullable=True)
    name = Column(String(1024), nullable=False)
    description = Column(Text, nullable=True)
    inchi = Column(Text, nullable=True)
    inchikey = Column(String(255), nullable=True)

    def __repr__(self):
        return self.name

    def as_bel(self) -> abundance:
        return abundance(namespace=MODULE_NAME, name=self.name, identifier=self.drugbank_id)

    def as_inchi_bel(self) -> abundance:
        # https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000383
        return abundance(namespace='inchi', name=self.inchi, identifier=self.inchi)

    def as_inchikey_bel(self) -> abundance:
        # https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000387
        return abundance(namespace='inchikey', name=self.inchikey, identifier=self.inchikey)

    def as_cas_bel(self) -> abundance:
        # https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000237
        return abundance(namespace='cas', name=self.cas_number, identifier=self.cas_number)


class DrugXref(Base):
    """Represents a drug's cross-reference to another database."""
    __tablename__ = DRUG_XREF_TABLE_NAME

    id = Column(Integer, primary_key=True)

    resource = Column(String(255), nullable=False)
    identifier = Column(String(255), nullable=False)

    drug_id = Column(Integer, ForeignKey(f'{DRUG_TABLE_NAME}.id'), nullable=False)
    drug = relationship(Drug, backref=backref('xrefs'))

    def __repr__(self):
        return f'{self.resource}:{self.identifier}'

    @classmethod
    def has_identifier(cls, resource, identifier):
        return and_(cls.resource == resource, cls.identifier == identifier)


class Patent(Base):
    """Represents a patent."""
    __tablename__ = PATENT_TABLE_NAME

    id = Column(Integer, primary_key=True)

    patent_id = Column(String(255), nullable=False, unique=True, index=True)
    country = Column(String(255), nullable=False)
    approved = Column(Date)
    expires = Column(Date)
    pediatric_extension = Column(Boolean)

    drugs = relationship(Drug, secondary=drug_patent, lazy='dynamic', backref=backref('patents', lazy='dynamic'))

    @staticmethod
    def filter_pk(country, patent_id):
        return and_(Patent.country == country, Patent.patent_id == patent_id)

    def to_json(self):
        return dict(
            patent_id=self.patent_id,
            country=self.country,
            approved_date=str(self.approved),
            expires_date=str(self.expires),
            pediatric_extension=self.pediatric_extension
        )

    @property
    def google_url(self) -> str:
        """Return the Google Patents URL of this patent."""
        return f'https://patents.google.com/patent/{PATENT_PREFIX_MAP[self.country]}{self.patent_id}'

    def download_pdfs(self, outdir: str):
        return download_google_patents(self.google_url, outdir)


class Alias(Base):
    """Represents an alias of a drug."""
    __tablename__ = ALIAS_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(Text)

    drug_id = Column(Integer, ForeignKey(f'{DRUG_TABLE_NAME}.id'), nullable=False)
    drug = relationship(Drug, backref=backref('aliases'))

    def __repr__(self):
        return self.name


class AtcCode(Base):
    """Represents an ATC code of a drug."""
    __tablename__ = ATC_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String(255))

    drug_id = Column(Integer, ForeignKey(f'{DRUG_TABLE_NAME}.id'), nullable=False)
    drug = relationship(Drug, backref=backref('atc_codes'))

    def __repr__(self):
        return self.name


class Group(Base):
    """Represents a drug group."""
    __tablename__ = GROUP_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String(255))

    drugs = relationship(Drug, secondary=drug_group, lazy='dynamic', backref=backref('groups', lazy='dynamic'))

    def __repr__(self):
        return self.name


class Category(Base):
    """Represents a drug category."""
    __tablename__ = CATEGORY_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String(255))
    mesh_id = Column(String(32))

    drugs = relationship(Drug, secondary=drug_category, lazy='dynamic', backref=backref('categories', lazy='dynamic'))

    def __repr__(self):
        return self.name


class Species(Base):
    """Represents a species."""
    __tablename__ = SPECIES_TABLE_NAME

    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True, index=True, nullable=False)

    def __repr__(self):
        return self.name


class Protein(Base):
    """Represents a protein."""
    __tablename__ = PROTEIN_TABLE_NAME

    id = Column(Integer, primary_key=True)

    species_id = Column(Integer, ForeignKey(f'{SPECIES_TABLE_NAME}.id'), nullable=False)
    species = relationship(Species)

    uniprot_id = Column(String(32))
    uniprot_accession = Column(String(32))
    name = Column(String(255))
    hgnc_id = Column(String(32))
    entrez_id = Column(String(32))

    def __repr__(self):
        return self.uniprot_id

    def as_bel_hgnc(self) -> protein:
        return protein(namespace='hgnc', identifier=self.hgnc_id)

    def as_bel(self) -> protein:
        """Serialize as a PyBEL node with the UniProt namespace."""
        return protein(namespace='uniprot', identifier=self.uniprot_id)


class Action(Base):
    """Represents the action a drug takes in a drug-protein interaction."""
    __tablename__ = ACTION_TABLE_NAME

    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True, index=True, nullable=False)


class Article(Base):
    """Represents an article in PubMed."""

    __tablename__ = ARTICLE_TABLE_NAME

    id = Column(Integer, primary_key=True)
    pubmed_id = Column(String(255), unique=True, nullable=False, index=True)


class DrugProteinInteraction(Base):
    """Represents an interaction between a drug and a protein."""
    __tablename__ = DRUG_PROTEIN_TABLE_NAME

    id = Column(Integer, primary_key=True)

    drug_id = Column(Integer, ForeignKey(f'{DRUG_TABLE_NAME}.id'))
    drug = relationship(Drug, backref=backref('protein_interactions'))

    protein_id = Column(Integer, ForeignKey(f'{PROTEIN_TABLE_NAME}.id'))
    protein = relationship(Protein, backref=backref('drug_interactions'))

    category = Column(String(32), nullable=False)  # target, enzyme, etc...
    known_action = Column(Boolean, nullable=False)

    actions = relationship(Action, secondary=dpi_action, lazy='dynamic',
                           backref=backref('drug_protein_interactions', lazy='dynamic'))

    articles = relationship(Article, secondary=dpi_article, lazy='dynamic',
                            backref=backref('drug_protein_interactions', lazy='dynamic'))

    def _add_to_graph(self, graph: BELGraph, u: BaseEntity, v: BaseEntity) -> Set[str]:
        """Return the set of keys used to add these edges."""
        return {
            graph.add_qualified_edge(
                u,
                v,
                relation=REGULATES,
                citation=article.pubmed_id,
                evidence='From DrugBank',
                annotations={
                    'bio2bel': MODULE_NAME,
                },
                object_modifier=activity(),
            )
            for article in self.articles
        }

    def add_to_graph(self, graph: BELGraph) -> Set[str]:
        """Add this interaction to the graph.

        :return: A set of the hashes of the edges that were added
        """
        # TODO update implementation to use actions

        drug_bel = self.drug.as_bel()
        protein_bel = self.protein.as_bel()

        return self._add_to_graph(graph, drug_bel, protein_bel)
