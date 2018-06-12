# -*- coding: utf-8 -*-

"""Suites for testing the populated database."""

from pybel import BELGraph
from tests.constants import PopulatedTemporaryCacheClassMixin

from pybel.constants import RELATION, CITATION, CITATION_TYPE, CITATION_REFERENCE, CITATION_TYPE_PUBMED
class TestPopulation(PopulatedTemporaryCacheClassMixin):
    """Tests the database is populated correctly."""

    def test_count(self):
        """Tests the correct number of drugs."""
        self.assertEqual(4, self.manager.count_drugs())
        self.assertLessEqual(23, self.manager.count_articles())

    def test_article(self):
        article = self.manager.get_article_by_pmid('10505536')
        self.assertIsNotNone(article)
        dpis = list(article.drug_protein_interactions)
        self.assertNotEqual(0, len(dpis))

    def test_bel(self):
        article = self.manager.get_article_by_pmid('10505536')

        drug_protein_interaction = article.drug_protein_interactions.all()[0]

        protein = drug_protein_interaction.protein
        self.assertEqual('P00734', protein.uniprot_id)

        drug = drug_protein_interaction.drug
        self.assertEqual('DB00001', drug.drugbank_id)
        self.assertEqual('Lepirudin', drug.name)

        graph = BELGraph()

        drug_protein_interaction.add_to_graph(graph)

        self.assertEqual(2, graph.number_of_nodes())
        self.assertEqual(1, graph.number_of_edges())

        _, _, data = list(graph.edges(data=True))[0]

        self.assertIn(CITATION, data)
        self.assertIn(CITATION_TYPE, data[CITATION])
        self.assertEqual(CITATION_TYPE_PUBMED, data[CITATION][CITATION_TYPE])
        self.assertIn(CITATION_REFERENCE, data[CITATION])
        self.assertEqual('10505536', data[CITATION][CITATION_REFERENCE])