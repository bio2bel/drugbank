# -*- coding: utf-8 -*-

"""Utilities for downloading patents from Google.

Code modified from original work by Alexander Esser.
"""

import os
import re
from typing import Optional, Set

import requests
from bs4 import BeautifulSoup

LINK_PATTERN = r"https?:\/\/patentimages\.storage\.googleapis\.com\/.+\/([A-z0-9]+\.pdf)"
LINK_RE = re.compile(LINK_PATTERN, re.IGNORECASE)

prefix_map = {
    'United States': 'US',
    'Canada': 'CA',
}


def download_google_patents(url: str, directory: str) -> Set[str]:
    """Crawls a list of URLs at patent.google.com and downloads the attached PDF documents

    :param url: The url (e.g., https://patents.google.com/patent/US5972916)
    :param directory: The output directory
    """
    rv = set()
    try:
        r = requests.get(url)
        data = r.text
        soup = BeautifulSoup(data, "html.parser")
        for link in soup.find_all("a"):
            target = link.get("href")
            link = _process_link(target, directory)
            if link:
                rv.add(link)

    except Exception as e:
        print("Could not download patent from {}: {}".format(url, str(e)))

    return rv


def _process_link(target, directory: str) -> Optional[str]:
    """Download the link if it fits the description and return it if it works."""
    m = LINK_RE.search(target)
    if not m:
        return

    outfile = os.path.join(directory, m.group(1))
    if os.path.exists(outfile):
        return target

    print(f"Downloading {target} to {outfile}")
    r2 = requests.get(target, stream=True)
    if r2.status_code != 200:
        return

    with open(outfile, 'wb') as f:
        for chunk in r2:
            f.write(chunk)

    return target
