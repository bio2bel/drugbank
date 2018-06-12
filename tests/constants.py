# -*- coding: utf-8 -*-

"""Constants for testing bio2bel_drugbank."""

import logging
import os

from bio2bel.testing import make_temporary_cache_class_mixin
from bio2bel_drugbank.manager import Manager

log = logging.getLogger(__name__)

dir_path = os.path.dirname(os.path.realpath(__file__))
resources_path = os.path.join(dir_path, 'resources')
test_xml_path = os.path.join(resources_path, 'test.xml')


class PopulatedTemporaryCacheClassMixin(make_temporary_cache_class_mixin(Manager)):
    """Create a test suite that has a populated database."""

    manager: Manager

    @classmethod
    def populate(cls):
        """Override the populate hook."""
        cls.manager.populate(url=test_xml_path)
