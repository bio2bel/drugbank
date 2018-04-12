# -*- coding: utf-8 -*-

import logging
import os

from bio2bel.testing import make_temporary_cache_class_mixin
from bio2bel_drugbank.manager import Manager

log = logging.getLogger(__name__)

dir_path = os.path.dirname(os.path.realpath(__file__))
resources_path = os.path.join(dir_path, 'resources')
test_xml_path = os.path.join(resources_path, 'test.xml')

TemporaryCacheClassMixin = make_temporary_cache_class_mixin(Manager)


class PopulatedTemporaryCacheClassMixin(TemporaryCacheClassMixin):
    @classmethod
    def populate(cls):
        cls.manager.populate(url=test_xml_path)
