# -*- coding: utf-8 -*-

"""Suites for testing the populated database."""

from tests.constants import PopulatedTemporaryCacheClassMixin


class TestPopulation(PopulatedTemporaryCacheClassMixin):
    """Tests the database is populated correctly."""

    def test_count(self):
        """Tests the correct number of drugs."""
        self.assertEqual(4, self.manager.count_drugs())
