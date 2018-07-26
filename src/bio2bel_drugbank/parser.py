# -*- coding: utf-8 -*-


import itertools as itt

import logging
import re
import time
from datetime import datetime
from xml.etree import ElementTree

from .constants import DRUGBANK_PATH

log = logging.getLogger(__name__)

ns = '{http://www.drugbank.ca}'
inchikey_template = f"{ns}calculated-properties/{ns}property[{ns}kind='InChIKey']/{ns}value"
inchi_template = f"{ns}calculated-properties/{ns}property[{ns}kind='InChI']/{ns}value"

pubmed_re = re.compile('pubmed/([0-9]+)')


def get_path(url=None):
    return url or DRUGBANK_PATH


def get_xml_root(url=None):
    """Get the XML parser root. Takes between 35-60 seconds.

    :param Optional[str] url: A custom URL for drugbank XML file
    :return:
    """
    url = get_path(url=url)

    log.info('parsing drugbank at %s', url)
    t = time.time()
    tree = ElementTree.parse(url)
    log.info('parsed drugbank in %.2f seconds', time.time() - t)

    return tree.getroot()


def extract_drug_info(drug_xml):
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
                'mesh_id': x.findtext(f'{ns}mesh-id')
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
                'identifier': x.findtext(f'{ns}identifier')
            }
            for x in drug_xml.findall(f"{ns}external-identifiers/{ns}external-identifier")
        ],
        'inchi': drug_xml.findtext(inchi_template),
        'inchikey': drug_xml.findtext(inchikey_template)
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

    for category, protein in _iterate_protein_stuff(drug_xml):
        protein_row = extract_protein_info(category, protein)
        if not protein_row:
            continue
        row['protein_interactions'].append(protein_row)

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
    row = {
        'category': category,
        'organism': protein.findtext(f'{ns}organism'),
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
        ]
    }

    # TODO generalize to all cross references?
    uniprot_ids = [polypep.text for polypep in protein.findall(
        f"{ns}polypeptide/{ns}external-identifiers/{ns}external-identifier[{ns}resource='UniProtKB']/{ns}identifier")]
    if len(uniprot_ids) != 1:
        return
    row['uniprot_id'] = uniprot_ids[0]

    hgnc_ids = [polypep.text for polypep in protein.findall(
        f"{ns}polypeptide/{ns}external-identifiers/{ns}external-identifier[{ns}resource='HUGO Gene Nomenclature Committee (HGNC)']/{ns}identifier")]
    if len(hgnc_ids) == 1:
        row['hgnc_id'] = hgnc_ids[0][len('HGNC:'):]

    return row
