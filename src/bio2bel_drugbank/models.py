# -*- coding: utf-8 -*-

"""Database models for bio2bel_drugbank"""

from sqlalchemy import Boolean, Column, Date, ForeignKey, Integer, String, Table, Text, and_
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import backref, relationship

__all__ = [
    'Base',
    'Type',
    'Drug',
    'Alias',
    'AtcCode',
    'Group',
    'Category',
    'Patent',
    'Xref',
]

Base = declarative_base()

TABLE_PREFIX = 'drugbank'
DRUG_TABLE_NAME = f'{TABLE_PREFIX}_drug'
TYPE_TABLE_NAME = f'{TABLE_PREFIX}_type'
ALIAS_TABLE_NAME = f'{TABLE_PREFIX}_alias'
ATC_TABLE_NAME = f'{TABLE_PREFIX}_atc'
CATEGORY_TABLE_NAME = f'{TABLE_PREFIX}_category'
GROUP_TABLE_NAME = f'{TABLE_PREFIX}_group'
DRUG_CATEGORY_TABLE_NAME = f'{TABLE_PREFIX}_drug_category'
DRUG_GROUP_TABLE_NAME = f'{TABLE_PREFIX}_drug_group'
PATENT_TABLE_NAME = f'{TABLE_PREFIX}_patent'
DRUG_PATENT_TABLE_NAME = f'{TABLE_PREFIX}_drug_patent'
XREF_TABLE_NAME = f'{TABLE_PREFIX}_xref'

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


class Type(Base):
    __tablename__ = TYPE_TABLE_NAME
    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True, index=True, nullable=False)

    def __repr__(self):
        return self.name


class Drug(Base):
    """Represents a chemical"""
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


class Xref(Base):
    """Represents a cross-reference to another database"""
    __tablename__ = XREF_TABLE_NAME

    id = Column(Integer, primary_key=True)

    resource = Column(String(255), nullable=False)
    identifier = Column(String(255), nullable=False)

    drug_id = Column(Integer, ForeignKey(f'{DRUG_TABLE_NAME}.id'), nullable=False)
    drug = relationship(Drug, backref=backref('xrefs'))

    def __repr__(self):
        return f'{self.resource}:{self.identifier}'


class Patent(Base):
    """Represents a patent"""
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


class Alias(Base):
    """Represents an alias of a drug"""
    __tablename__ = ALIAS_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(Text)

    drug_id = Column(Integer, ForeignKey(f'{DRUG_TABLE_NAME}.id'), nullable=False)
    drug = relationship(Drug, backref=backref('aliases'))

    def __repr__(self):
        return self.name


class AtcCode(Base):
    """Represents an ATC code of a drug"""
    __tablename__ = ATC_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String(255))

    drug_id = Column(Integer, ForeignKey(f'{DRUG_TABLE_NAME}.id'), nullable=False)
    drug = relationship(Drug, backref=backref('atc_codes'))

    def __repr__(self):
        return self.name


class Group(Base):
    __tablename__ = GROUP_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String(255))

    drugs = relationship(Drug, secondary=drug_group, lazy='dynamic', backref=backref('groups', lazy='dynamic'))

    def __repr__(self):
        return self.name


class Category(Base):
    __tablename__ = CATEGORY_TABLE_NAME

    id = Column(Integer, primary_key=True)

    name = Column(String(255))
    mesh_id = Column(String(32))

    drugs = relationship(Drug, secondary=drug_category, lazy='dynamic', backref=backref('categories', lazy='dynamic'))

    def __repr__(self):
        return self.name
