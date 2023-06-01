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
from src.eda import ExploratoryAnalyzer


class TestDataCleanerCollectData(unittest.TestCase):
    ''' Tests the DataCleaner 'collect_data()' function. '''
    
    def setUp(self):
        self.cleaner = DataCleaner(input_path=Path('mock_data_for_tests/data_cleaner/collect_data_function'))

    def test_data_collection_raises_error(self):
        ''' Empty file, function should raise an error. '''

        print("\nTesting that collecting from a directory with no files returns error.")
        # Change path to test for empty directory. 
        self.cleaner.input_path = Path('mock_data_for_tests/data_cleaner/empty_file')
        self.assertRaises(FileNotFoundError, self.cleaner.collect_data)
    
    def test_data_collection_with_column_errors(self):
        ''' Data frame does not have the correct info, function should remove it. '''

        print("\nTesting that collecting from a directory with column errors returns error.")
        collect_data = self.cleaner.collect_data()
        self.assertEqual(len(collect_data), 2)

    def test_data_collection_reindexes_correctly(self):
        ''' Columns are out of order, function should reorder the columns. '''

        print("\nTesting extra columns are removed and are indexed into the correct order.")
        collect_data = self.cleaner.collect_data()
        self.assertEqual(len(collect_data['bad_column_order.tsv'].columns), len(self.cleaner.column_order))

    
class TestDataCleanerRemoveBad(unittest.TestCase):
    ''' Tests the DataCleaner 'remove_bad_reads()' function. '''

    def setUp(self) -> None:
        self.cleaner = DataCleaner(input_path=Path('mock_data_for_tests/data_cleaner/remove_bad_reads_function'))

    def test_bad_sequences_are_removed(self):
        ''' File contains bad sequences, function should remove them. '''

        print("\nTesting that 'remove_bad_reads()' removes sequence errors from the data frame.")
        collected_data = self.cleaner.collect_data()
        bad_seqs_removed = self.cleaner.remove_bad_reads(collected_data)
        self.assertEqual(bad_seqs_removed['contains_bad_reads.tsv']['AASeq'].values.size, 2)


class TestDataCleanerAssertSize(unittest.TestCase):
    ''' Tests the DataCleaner 'assert_size()' function. '''

    def setUp(self):
        self.cleaner = DataCleaner(input_path=Path('mock_data_for_tests/data_cleaner/assert_size_function'))

    def test_bad_sizes_are_removed(self):
        ''' File contains sequences that are too large and too small, function should remove them. '''
    
        print("\nTesting that 'assert_size()' removes non-biologically normal seqs from the data frame.")
        collected_data = self.cleaner.collect_data()
        bad_seqs_removed = self.cleaner.remove_bad_reads(collected_data)
        abnormal_seqs_removed = self.cleaner.assert_size(bad_seqs_removed)
        self.assertEqual(abnormal_seqs_removed['contains_abnormal_reads.tsv']['AASeq'].values.size, 2)


class TestDataCleanerJoinFiles(unittest.TestCase):
    ''' Tests DataCleaner class 'join_files()' function. '''

    def setUp(self):
        self.cleaner = DataCleaner(input_path=Path('mock_data_for_tests/data_cleaner/join_files_function'))

    def test_not_supported_disease_type(self):
        ''' File name is not part of the training set, function should raise and error. '''
        
        print("\nTesting that an error is raised when the file name is not part of the training set.")
        collected_data = self.cleaner.collect_data()
        bad_seqs_removed = self.cleaner.remove_bad_reads(collected_data)
        abnormal_seqs_removed = self.cleaner.assert_size(bad_seqs_removed)
        self.assertRaises(DiseaseNotSupportedError, self.cleaner.join_files, abnormal_seqs_removed)

class TestExploratoryAnalysis(unittest.TestCase):
    ''' Tests the EDA Class. '''
    
    def setUp(self):
        self.eda = ExploratoryAnalyzer(collection_path=Path('mock_data_for_tests/eda'))
    
    def test_count_seqs_per_disease(self):
        ''' Tests the 'count_seqs_per_disease()' function. '''

        print("\nTesting that the 'count_seqs_per_disease()' function returns the correct number of sequences per disease.")
        collected_data = self.eda.collect_files()
        self.eda.count_sequences_per_disease(collected_data)
        self.assertEqual(self.eda.sequence_counts['test_1'], 3)
        self.assertEqual(self.eda.sequence_counts['test_2'], 5)
        self.assertEqual(self.eda.sequence_counts['test_3'], 2)

    def test_sequence_length_metrics(self):
        print("\nTesting that the 'sequence_length_metrics()' function returns the correct sequence length metrics.")
        collected_data = self.eda.collect_files()
        self.eda.count_sequences_per_disease(collected_data)
        self.eda.count_sequence_length_metrics(collected_data)
        self.assertEqual(list(self.eda.seq_len_metrics['test_1']['len_counts']), [1.0])
        self.assertEqual(list(self.eda.seq_len_metrics['test_2']['len_counts']), [0.4, 0.2, 0.2, 0.2])
        self.assertEqual(list(self.eda.seq_len_metrics['test_3']['len_counts']), [0.5, 0.5])

    def test_gene_usage(self):
        print("\nTesting that the 'count_gene_usage()' function returns the correct gene usage.")

        collected_data = self.eda.collect_files()
        self.eda.count_gene_usage(collected_data)
        # Test 1.
        self.assertEqual(list(self.eda.gene_usage['Vregion']['test_1']), [0.6666666666666666, 0.3333333333333333])
        self.assertEqual(list(self.eda.gene_usage['Dregion']['test_1']), [0.6666666666666666, 0.3333333333333333])
        self.assertEqual(list(self.eda.gene_usage['Jregion']['test_1']), [1.0])
        # Test 2.
        self.assertEqual(list(self.eda.gene_usage['Vregion']['test_2']), [0.2, 0.2, 0.2, 0.2, 0.2])
        self.assertEqual(list(self.eda.gene_usage['Dregion']['test_2']), [1.0])
        self.assertEqual(list(self.eda.gene_usage['Jregion']['test_2']), [0.4, 0.2, 0.2, 0.2])
        # Test 3.
        self.assertEqual(list(self.eda.gene_usage['Vregion']['test_3']), [0.5, 0.5])
        self.assertEqual(list(self.eda.gene_usage['Dregion']['test_3']), [1.0])
        self.assertEqual(list(self.eda.gene_usage['Jregion']['test_3']), [0.5, 0.5])

    def test_jaccard_index(self):
        print("\nTesting that the 'jaccard_index()' function returns the correct Jaccard Index.")
        collected_data = self.eda.collect_files()
        self.eda.compute_jaccard_index(collected_data)
        self.assertEqual(self.eda.jaccard_index['test_1'],  {'test_2': 0.0, 'test_3': 0.0})
        self.assertEqual(self.eda.jaccard_index['test_2'], {'test_1': 0.0, 'test_3': 0.16666666666666666})
        self.assertEqual(self.eda.jaccard_index['test_3'], {'test_1': 0.0, 'test_2': 0.16666666666666666})

    def test_aa_occurrences(self):
        print("\nTesting that the 'count_aa_occurrences()' function returns the correct AA occurrences.")
        collected_data = self.eda.collect_files()
        self.eda.count_aa_occurrences(collected_data)
        print(list(self.eda.aa_counts['test_1'][1]))
        # Test 1.
        self.assertEqual(list(self.eda.aa_counts['test_1'][1]), [1.0])
        self.assertEqual(list(self.eda.aa_counts['test_1'][5]), [0.6666666666666666, 0.3333333333333333])
        self.assertEqual(list(self.eda.aa_counts['test_1'][9]), [0.6666666666666666, 0.3333333333333333])
        # Test 2.
        self.assertEqual(list(self.eda.aa_counts['test_2'][1]), [1.0])
        self.assertEqual(list(self.eda.aa_counts['test_2'][5]), [0.4, 0.2, 0.2, 0.2])
        self.assertEqual(list(self.eda.aa_counts['test_2'][9]), [0.6, 0.2, 0.2])
        # Test 3.
        self.assertEqual(list(self.eda.aa_counts['test_3'][1]), [1.0])
        self.assertEqual(list(self.eda.aa_counts['test_3'][5]), [0.5, 0.5])
        self.assertEqual(list(self.eda.aa_counts['test_3'][9]), [0.5, 0.5])
        


# class PreprocessingTester(unittest.TestCase):
#     ''' Tests the Preprocessing Class. '''
#     pass



unittest.main()