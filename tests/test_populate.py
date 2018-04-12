# -*- coding: utf-8 -*-

from tests.constants import PopulatedTemporaryCacheClassMixin


class TestPopulation(PopulatedTemporaryCacheClassMixin):

    def test_count(self):
        self.assertEqual(4, self.manager.count_drugs())
