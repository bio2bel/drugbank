# -*- coding: utf-8 -*-

"""Utilities for downloading patents from Google.

Code modified from original work by Alexander Esser.
"""

import os
import re

import requests
from bs4 import BeautifulSoup

LINK_PATTERN = "https?:\/\/patentimages\.storage\.googleapis\.com\/.+\/([A-z0-9]+\.pdf)"
LINK_RE = re.compile(LINK_PATTERN, re.IGNORECASE)

prefix_map = {
    'United States': 'US',
    'Canada': 'CA',
}


def download_google_patents(url: str, directory:str):
    """Crawls a list of URLs at patent.google.com and downloads the attached PDF documents

    :param url: The url (e.g., https://patents.google.com/patent/US5972916)
    :param directory: The output directory
    """
    try:
        r = requests.get(url)
        data = r.text
        soup = BeautifulSoup(data, "html.parser")
        for link in soup.find_all("a"):
            _process_link(link, directory)

    except Exception as e:
        print("Could not download patent from {}: {}".format(url, str(e)))


def _process_link(link, directory:str):
    target = link.get("href")
    m = LINK_RE.search(target)
    if not m:
        return

    r2 = requests.get(target, stream=True)
    if r2.status_code != 200:
        return

    outfile = os.path.join(directory, m.group(1))
    print("Downloading {} to {}".format(target, outfile))
    with open(outfile, 'wb') as f:
        for chunk in r2:
            f.write(chunk)
