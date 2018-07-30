# -*- coding: utf-8 -*-

import os

from bio2bel.utils import get_connection, get_data_dir

VERSION = '0.1.0'

MODULE_NAME = 'drugbank'
DATA_DIR = get_data_dir(MODULE_NAME)
DEFAULT_CACHE_CONNECTION = get_connection(MODULE_NAME)

DRUGBANK_URL = 'https://www.drugbank.ca/releases/5-1-0/downloads/all-full-database'
DRUGBANK_PATH = os.path.join(DATA_DIR, 'full database.xml')
