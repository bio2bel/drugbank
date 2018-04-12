# -*- coding: utf-8 -*-


import itertools as itt
import logging
import re
import time
from datetime import datetime
from xml.etree import ElementTree

from tqdm import tqdm

from .constants import DRUGBANK_PATH

log = logging.getLogger(__name__)

ns = '{http://www.drugbank.ca}'
inchikey_template = f"{ns}calculated-properties/{ns}property[{ns}kind='InChIKey']/{ns}value"
inchi_template = f"{ns}calculated-properties/{ns}property[{ns}kind='InChI']/{ns}value"


def get_path():
    return DRUGBANK_PATH


def get_xml_root(url=None):
    """Get the XML parser root

    Took about 35 seconds

    :param Optional[str] url: A custom URL for drugbank XML file
    :return:
    """
    log.info('parsing drugbank')
    t = time.time()
    tree = ElementTree.parse(url or DRUGBANK_PATH)
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
                'resource':  x.findtext(f'{ns}resource'),
                'identifier':  x.findtext(f'{ns}identifier')
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

    return row


_categories = ['target', 'enzyme', 'carrier', 'transporter']


def _iterate_protein_stuff(drug):
    drugbank_id = drug.findtext(ns + "drugbank-id[@primary='true']")
    for category in _categories:
        proteins = drug.findall(f'{ns}{category}s/{ns}{category}')
        for protein in proteins:
            yield drugbank_id, category, protein


def extract_protein_info(drugbank_id, category, protein):
    row = {
        'drugbank_id': drugbank_id,
        'category': category,
        'organism': protein.findtext(f'{ns}organism'),
        'known_action': protein.findtext(f'{ns}known-action')
    }

    actions = protein.findall(f'{ns}actions/{ns}action')
    row['actions'] = '|'.join(action.text for action in actions)

    uniprot_ids = [polypep.text for polypep in protein.findall(
        f"{ns}polypeptide/{ns}external-identifiers/{ns}external-identifier[{ns}resource='UniProtKB']/{ns}identifier")]
    if len(uniprot_ids) != 1:
        return
    row['uniprot_id'] = uniprot_ids[0]

    ref_text = protein.findtext(f"{ns}references[@format='textile']")
    pmids = re.findall(r'pubmed/([0-9]+)', ref_text)
    row['pubmed_ids'] = '|'.join(pmids)

    return row


def _iterate_protein_info(root):
    for drug in root:
        for drugbank_id, category, protein in _iterate_protein_stuff(drug):
            row = extract_protein_info(drugbank_id, category, protein)
            if row:
                yield row


def get_drug_info(root):
    log.info('getting all drug info')
    return [
        extract_drug_info(drug)
        for drug in tqdm(root)
    ]
