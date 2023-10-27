import unittest
import pandas as pd
from parse_BLAST_results import filter_hits
#import scripts.parse_BLAST_results


class test_filter_hits(unittest.TestCase):
    def setUp(self):
        self.input_file = 'tests/inputs/sample_22_viral_contigs_blasted.txt'

    def test_filter_hits(self):
        self.assertIsInstance(filter_hits(self.input_file, 'Viruses'), pd.DataFrame, 'filter_hits() returns a DataFrame')

unittest.main()