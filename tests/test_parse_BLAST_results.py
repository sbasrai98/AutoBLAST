import unittest
import pandas as pd
from scripts.parse_BLAST_results import *

class test_filter_hits(unittest.TestCase):
    def setUp(self):
        self.input_file = 'tests/inputs/sample_22_viral_contigs_blasted.txt'

    def test_filter_hits(self):
        parsed_hits = filter_hits(self.input_file, 'Viruses')
        self.assertIsInstance(parsed_hits, pd.DataFrame, 'filter_hits() returns a DataFrame')
        self.assertEqual(parsed_hits.shape[0], 70, 'parsed_hits DataFrame has correct number of rows')
        self.assertEqual(parsed_hits.shape[1], 8, 'parsed_hits DataFrame has correct number of columns')
        