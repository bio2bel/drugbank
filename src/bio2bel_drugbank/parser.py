# -*- coding: utf-8 -*-

"""Parsers for Bio2BEL DrugBank."""

import itertools as itt
import logging
import re
from datetime import datetime
from typing import Mapping, Optional
from xml.etree import ElementTree

from drugbank_downloader import parse_drugbank
from tqdm import tqdm

import pyobo.config


log = logging.getLogger(__name__)

ns = '{http://www.drugbank.ca}'
inchikey_template = f"{ns}calculated-properties/{ns}property[{ns}kind='InChIKey']/{ns}value"
inchi_template = f"{ns}calculated-properties/{ns}property[{ns}kind='InChI']/{ns}value"
smiles_template = f"{ns}calculated-properties/{ns}property[{ns}kind='SMILES']/{ns}value"

pubmed_re = re.compile('pubmed/([0-9]+)')


def get_xml_root(path: Optional[str] = None):
    """Get the XML parser root.

    Takes between 35-60 seconds.

    :param path: A custom URL for DrugBank XML file
    """
    if path:
        return ElementTree.parse(path).getroot()
    return parse_drugbank(
        username=pyobo.config.get_config('drugbank_username'),
        password=pyobo.config.get_config('drugbank_password'),
    )


def extract_drug_info(drug_xml: ElementTree.Element):
    """Extract information from an XML element representing a drug."""
    assert drug_xml.tag == f'{ns}drug'

    row = {
        'type': drug_xml.get('type'),
        'drugbank_id': drug_xml.findtext(f"{ns}drugbank-id[@primary='true']"),
        'cas_number': drug_xml.findtext(f"{ns}cas-number"),
        'name': drug_xml.findtext(f"{ns}name"),
        'description': drug_xml.findtext(f"{ns}description"),
        'groups': [
            group.text
            for group in drug_xml.findall(f"{ns}groups/{ns}group")
        ],
        'atc_codes': [
            code.get('code')
            for code in drug_xml.findall(f"{ns}atc-codes/{ns}atc-code")
        ],
        'categories': [
            {
                'name': x.findtext(f'{ns}category'),
                'mesh_id': x.findtext(f'{ns}mesh-id'),
            }
            for x in drug_xml.findall(f"{ns}categories/{ns}category")
        ],
        'patents': [
            {
                'patent_id': x.findtext(f'{ns}number'),
                'country': x.findtext(f'{ns}country'),
                'approved': datetime.strptime(x.findtext(f'{ns}approved'), '%Y-%m-%d'),
                'expires': datetime.strptime(x.findtext(f'{ns}expires'), '%Y-%m-%d'),
                'pediatric_extension': x.findtext(f'{ns}pediatric-extension') != 'false',

            }
            for x in drug_xml.findall(f"{ns}patents/{ns}patent")
        ],
        'xrefs': [
            {
                'resource': x.findtext(f'{ns}resource'),
                'identifier': x.findtext(f'{ns}identifier'),
            }
            for x in drug_xml.findall(f"{ns}external-identifiers/{ns}external-identifier")
        ],
        'inchi': drug_xml.findtext(inchi_template),
        'inchikey': drug_xml.findtext(inchikey_template),
        'smiles': drug_xml.findtext(smiles_template),
    }

    # Add drug aliases
    aliases = {
        elem.text.strip() for elem in
        itt.chain(
            drug_xml.findall(f"{ns}international-brands/{ns}international-brand"),
            drug_xml.findall(f"{ns}synonyms/{ns}synonym[@language='English']"),
            drug_xml.findall(f"{ns}international-brands/{ns}international-brand"),
            drug_xml.findall(f"{ns}products/{ns}product/{ns}name"),
        )
        if elem.text.strip()
    }
    aliases.add(row['name'])
    row['aliases'] = aliases

    row['protein_interactions'] = []
    row['protein_group_interactions'] = []

    for category, protein in _iterate_protein_stuff(drug_xml):
        target_row = extract_protein_info(category, protein)
        if not target_row:
            continue
        row['protein_interactions'].append(target_row)

    return row


_categories = ['target', 'enzyme', 'carrier', 'transporter']


def _iterate_protein_stuff(drug_xml):
    """

    :param drug_xml:
    :return: iterates over pairs of category and protein xml
    """
    for category in _categories:
        proteins = drug_xml.findall(f'{ns}{category}s/{ns}{category}')
        for protein in proteins:
            yield category, protein


def extract_protein_info(category, protein):
    # FIXME Differentiate between proteins and protein groups/complexes
    polypeptides = protein.findall(f'{ns}polypeptide')

    if len(polypeptides) == 0:
        protein_type = 'none'
    elif len(polypeptides) == 1:
        protein_type = 'single'
    else:
        protein_type = 'group'

    row = {
        'target': {
            'type': protein_type,
            'category': category,
            'known_action': protein.findtext(f'{ns}known-action'),
            'name': protein.findtext(f'{ns}name'),
            'actions': [
                action.text
                for action in protein.findall(f'{ns}actions/{ns}action')
            ],
            'articles': [
                pubmed_element.text
                for pubmed_element in protein.findall(f'{ns}references/{ns}articles/{ns}article/{ns}pubmed-id')
                if pubmed_element.text
            ],
        },
        'polypeptides': [
            {
                'name': polypeptide.findtext(f'{ns}name'),
                'hgnc_symbol': polypeptide.findtext(f'{ns}gene-name'),
                'hgnc_id': polypeptide.findtext(
                    f"{ns}external-identifiers/{ns}external-identifier[{ns}resource='HUGO Gene Nomenclature Committee (HGNC)']/{ns}identifier",
                )[len('HGNC:'):],
                'uniprot_id': polypeptide.findtext(
                    f"{ns}external-identifiers/{ns}external-identifier[{ns}resource='UniProtKB']/{ns}identifier",
                ),
                'uniprot_accession': polypeptide.findtext(
                    f"{ns}external-identifiers/{ns}external-identifier[{ns}resource='UniProt Accession']/{ns}identifier",
                ),
                'organism': polypeptide.findtext(f'{ns}organism'),
                'taxonomy': polypeptide.find(f'{ns}organism').attrib['ncbi-taxonomy-id'],
            }
            for polypeptide in polypeptides
        ],
    }

    return row


def get_pubchem_to_drugbank(path=None) -> Mapping[str, str]:
    """Get a mapping from PubChem Substances to DrugBank identifiers."""
    rv = {}
    root = get_xml_root(path=path)
    for drug_xml in tqdm(root, desc='Drugs'):
        drug = extract_drug_info(drug_xml)
        drugbank_id = drug['drugbank_id']
        for xref in drug['xrefs']:
            if xref['resource'] == 'PubChem Substance':
                rv[xref['identifier']] = drugbank_id
                break
        else:
            print(f'could not find pubchem for {drugbank_id}')

    return rv


def main():
    import json
    import os

    x = get_pubchem_to_drugbank()
    with open(os.path.expanduser('~/Desktop/pubchem_to_drugbank.json'), 'w') as f:
        json.dump(x, f, indent=2)


if __name__ == '__main__':
    main()
