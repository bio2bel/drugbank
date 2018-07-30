# -*- coding: utf-8 -*-

"""
Installation
------------
``bio2bel_drugbank`` can be installed easily from `PyPI <https://pypi.python.org/pypi/bio2bel_drugbank>`_ with
the following code in your favorite terminal:

.. code-block:: sh

    $ python3 -m pip install bio2bel_drugbank

or from the latest code on `GitHub <https://github.com/bio2bel/drugbank>`_ with:

.. code-block:: sh

    $ python3 -m pip install git+https://github.com/bio2bel/drugbank.git@master

Setup
-----
.. warning:: DrugBank requires a bit of downloading and file organization. Will be documented soon.

Python REPL
~~~~~~~~~~~
.. code-block:: python

    >>> import bio2bel_drugbank
    >>> drugbank_manager = bio2bel_drugbank.Manager()
    >>> drugbank_manager.populate()

Command Line Utility
~~~~~~~~~~~~~~~~~~~~
.. code-block:: bash

    bio2bel_drugbank populate
"""

from . import manager, models, utils
from .manager import *
from .models import *
from .utils import *

__all__ = manager.__all__ + models.__all__ + utils.__all__

__version__ = '0.1.0-dev'

__title__ = 'bio2bel_drugbank'
__description__ = "A package for converting DrugBank to BEL"
__url__ = 'https://github.com/bio2bel/drugbank'

__author__ = 'Charles Tapley Hoyt'
__email__ = 'charles.hoyt@scai.fraunhofer.de'

__license__ = 'MIT License'
__copyright__ = 'Copyright (c) 2018 Charles Tapley Hoyt'
