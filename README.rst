Bio2BEL DrugBank |build| |coverage| |documentation|
===================================================
Converts DrugBank 5.0 to BEL.

Installation |pypi_version| |python_versions| |pypi_license|
------------------------------------------------------------
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

Citations
---------
- Wishart, D. S., *et al.* (2018). `DrugBank 5.0: a major update to the DrugBank database for 2018
  <https://doi.org/10.1093/nar/gkx1037>`_. Nucleic Acids Research, 46(D1), D1074–D1082.

.. |build| image:: https://travis-ci.org/bio2bel/drugbank.svg?branch=master
    :target: https://travis-ci.org/bio2bel/drugbank
    :alt: Build Status

.. |documentation| image:: http://readthedocs.org/projects/bio2bel-drugbank/badge/?version=latest
    :target: http://bio2bel.readthedocs.io/projects/drugbank/en/latest/?badge=latest
    :alt: Documentation Status

.. |pypi_version| image:: https://img.shields.io/pypi/v/bio2bel_drugbank.svg
    :alt: Current version on PyPI

.. |coverage| image:: https://codecov.io/gh/bio2bel/drugbank/coverage.svg?branch=master
    :target: https://codecov.io/gh/bio2bel/drugbank?branch=master
    :alt: Coverage Status

.. |climate| image:: https://codeclimate.com/github/bio2bel/drugbank/badges/gpa.svg
    :target: https://codeclimate.com/github/bio2bel/drugbank
    :alt: Code Climate

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/bio2bel_drugbank.svg
    :alt: Stable Supported Python Versions

.. |pypi_license| image:: https://img.shields.io/pypi/l/bio2bel_drugbank.svg
    :alt: MIT License
