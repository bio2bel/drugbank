# -*- coding: utf-8 -*-

from . import manager, models
from .manager import *
from .models import *

__all__ = (manager.__all__ + models.__all__)

__version__ = '0.0.1-dev'

__title__ = 'bio2bel_drugbank'
__description__ = "A package for converting DrugBank to BEL"
__url__ = 'https://github.com/bio2bel/drugbank'

__author__ = 'Charles Tapley Hoyt'
__email__ = 'charles.hoyt@scai.fraunhofer.de'

__license__ = 'MIT License'
__copyright__ = 'Copyright (c) 2018 Charles Tapley Hoyt'
