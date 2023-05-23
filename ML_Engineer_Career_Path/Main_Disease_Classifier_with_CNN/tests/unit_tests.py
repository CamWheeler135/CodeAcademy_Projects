''' File contains the unit tests. '''

# Module Imports.
import unittest
from pathlib import Path
import pandas as pd

# Allows the tests to get to the SRC file. 
import sys
sys.path.append('../')

# # Local Imports.
from src.data_cleaning import DataCleaner
from src.errors import FileNotFoundError, DiseaseNotSupportedError


class DataCleanerCollectData(unittest.TestCase):
    ''' Tests the DataCleaner 'collect_data()' function. '''
    
    def setUp(self):
        self.cleaner = DataCleaner(input_path=Path('mock_data_for_tests/data_cleaner/collect_data_function'))

    def test_data_collection_raises_error(self):
        ''' Empty file, function should raise an error. '''

        print("\nTesting that collecting from an empty directory returns error.")
        # Change path to test for empty directory. 
        self.cleaner.input_path = Path('mock_data_for_tests/data_cleaner/empty_file')
        self.assertRaises(FileNotFoundError, self.cleaner.collect_data)
    
    def test_data_collection_with_column_errors(self):
        ''' Data frame does not have the correct info, function should remove it. '''

        print("\nTesting that data frames without the needed information are removed")
        collect_data = self.cleaner.collect_data()
        self.assertEqual(len(collect_data), 2)

    def test_data_collection_reindexes_correctly(self):
        ''' Columns are out of order, function should reorder the columns. '''

        print("\nTesting extra columns are removed and are indexed into the correct order.")
        collect_data = self.cleaner.collect_data()
        self.assertEqual(len(collect_data['bad_column_order.tsv'].columns), len(self.cleaner.column_order))

    
class DataCleanerRemoveBad(unittest.TestCase):
    ''' Tests the DataCleaner 'remove_bad_reads()' function. '''

    def setUp(self) -> None:
        self.cleaner = DataCleaner(input_path=Path('mock_data_for_tests/data_cleaner/remove_bad_reads_function'))

    def test_bad_sequences_are_removed(self):
        ''' File contains bad sequences, function should remove them. '''

        print("\nTesting that bad sequences reads are removed from the data frame. ")
        collected_data = self.cleaner.collect_data()
        bad_seqs_removed = self.cleaner.remove_bad_reads(collected_data)
        self.assertEqual(bad_seqs_removed['contains_bad_reads.tsv']['AASeq'].values.size, 2)


class DataCleanerAssertSize(unittest.TestCase):
    ''' Tests the DataCleaner 'assert_size()' function. '''

    def setUp(self):
        self.cleaner = DataCleaner(input_path=Path('mock_data_for_tests/data_cleaner/assert_size_function'))

    def test_bad_sizes_are_removed(self):
        ''' File contains sequences that are too large and too small, function should remove them. '''
    
        print("\nTesting that abnormal sequences are removed from the data frame. ")
        collected_data = self.cleaner.collect_data()
        bad_seqs_removed = self.cleaner.remove_bad_reads(collected_data)
        abnormal_seqs_removed = self.cleaner.assert_size(bad_seqs_removed)
        self.assertEqual(abnormal_seqs_removed['contains_abnormal_reads.tsv']['AASeq'].values.size, 2)


class DataCleanerJoinFiles(unittest.TestCase):
    ''' Tests DataCleaner class 'join_files()' function. '''

    def setUp(self):
        self.cleaner = DataCleaner(input_path=Path('mock_data_for_tests/data_cleaner/join_files_function'))

    def test_not_supported_disease_type(self):
        ''' File name is not part of the training set, function should raise and error. '''
        
        print("\nTesting that currently not support disease types cause an exception. ")
        collected_data = self.cleaner.collect_data()
        bad_seqs_removed = self.cleaner.remove_bad_reads(collected_data)
        abnormal_seqs_removed = self.cleaner.assert_size(bad_seqs_removed)
        self.assertRaises(DiseaseNotSupportedError, self.cleaner.join_files, abnormal_seqs_removed)

# class ExploratoryAnalysisTester(unittest.TestCase):
#     ''' Tests the EDA Class. '''
#     pass


# class PreprocessingTester(unittest.TestCase):
#     ''' Tests the Preprocessing Class. '''
#     pass



unittest.main()