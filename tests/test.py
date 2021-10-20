"""imports"""
import unittest
from whispr import *

class MyTests(unittest.TestCase):
    """unit tests for utility functions"""

    def test_generateVolumeTable(self):

        #plate_type, vol_table, source_plate_layout, output_layout

        #CSV file with source plate info. Columns are label, well, concentration, volume
        source_plate = 'source_plates/ECHO_source_plate.xlsx'
        source_plate_df = pd.read_excel(source_plate, index_col = 0) 

        #CSV file with each reaction as rows and the volume of each input (columns) added to reach final reaction volume
        mixing_table = 'mixing_tables/mixing_table_conc_test.csv' 
        mixing_table_df = pd.read_csv(mixing_table, index_col = 0, dtype = str)

        vol_table = generateVolumeTable(mixing_table_df, source_plate_df) # TODO : real test
        self.assertTrue(True)
