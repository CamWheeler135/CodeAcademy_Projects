''' File contains the unit tests. '''

# Module Imports.
import unittest
from pathlib import Path
import pandas as pd
import os

# Allows the tests to get to the SRC file. 
import sys
sys.path.append('../')

# Local Imports.
from src.data_cleaning import DataHandler, DataCleaner
from src.errors import FileNotFoundError, DiseaseNotSupportedError
from src.eda import ExploratoryAnalyzer
from src.preprocessing import KmerPreprocessor

class TestDataHandler(unittest.TestCase):
    ''' 
    Tests the DataHandler class.
    - Instances are not made using the setUp function as we need to different paths for each of the tests.
    '''

    def test_load_raw_tcrdb(self):
        ''' Tests that the handler will load the data into dataframes correctly. '''
        self.raw_handler = DataHandler(collection_path=Path('mock_data_for_tests/data_handler/raw_files'),
                                       save_path=Path('This_is_a_mock_path')) # Add a save path to avoid 'warning' message.
        repertoire_dict = self.raw_handler.collect_raw_tcrdb_files()
        self.assertEqual(len(repertoire_dict), 2)

    def test_no_file_found(self):
        ''' Tests that the handler will raise a FileNotFoundError if the file is not found. '''
        self.empty_handler = DataHandler(collection_path=Path('mock_data_for_tests/data_handler/empty_folder'),
                                         save_path=Path('This_is_a_mock_path')) # Add a save path to avoid 'warning' message.
        self.assertRaises(FileNotFoundError, self.empty_handler.collect_raw_tcrdb_files)
   
    def test_project_files_are_handled(self):
        ''' Tests that the handler will load the project and samples files correctly.'''
        self.project_file_handler = DataHandler(collection_path=Path('mock_data_for_tests/data_handler/raw_files'),
                                                save_path=Path('This_is_a_mock_path')) # Add a save path to avoid 'warning' message.
        repertoire_dict = self.project_file_handler.collect_raw_tcrdb_files()
        project_file_df_columns = repertoire_dict['PRJ12345.tsv'].columns
        sample_file_df_columns = repertoire_dict['sample12345.tsv'].columns
        self.assertEqual(len(project_file_df_columns), 4)
        self.assertEqual(len(sample_file_df_columns), 4)

    '''TODO'''
    # We need to catch the error where not all the information is in the TSV file.

    def test_load_csv(self):
        ''' Tests that the handler will load the clean csv files correctly. '''
        self.csv_handler = DataHandler(collection_path=Path('mock_data_for_tests/data_handler/csv_files'),
                                       save_path=Path('This_is_a_mock_path')) # Add a save path to avoid 'warning' message.
        csv_dict = self.csv_handler.collect_csv_files()
        file_1_columns = csv_dict['csv1'].columns
        file_2_columns = csv_dict['csv2'].columns
        self.assertEqual(len(csv_dict), 2) # Check that the files are all loaded.
        self.assertEqual(len(file_1_columns), 4)
        self.assertEqual(len(file_2_columns), 4)
    
    def test_save_csv(self):
        ''' Tests that the handler will save the csv files correctly. '''
        self.save_csv_handler = DataHandler(collection_path=Path('mock_data_for_tests/data_handler/csv_files'),
                                       save_path=Path('mock_data_for_tests/data_handler/save_csv_test/saved_csv')) # Normal construction of the Handler.
        csv_dict = self.save_csv_handler.collect_csv_files()
        file_for_df = csv_dict['csv1']
        self.save_csv_handler.save_as_csv(file_for_df)
        files_in_dir = os.listdir(os.path.split(self.save_csv_handler.save_path)[0])
        self.assertEqual(len(files_in_dir), 1)


class TestDataCleaner(unittest.TestCase):
    ''' Tests the DataCleaner class. '''

    def setUp(self):
        handler = DataHandler(collection_path=Path('mock_data_for_tests/data_cleaner'),
                                   save_path=Path('This_is_a_mock_path')) # Add a save path to avoid 'warning' message
        self.cleaner = DataCleaner(handler=handler)

    def test_remove_bad_reads_(self):
        ''' Tests that sequences with errors are removed. '''
        repertoire_dict = self.cleaner.collect_files()
        removed_bad_reads = self.cleaner.remove_bad_reads(repertoire_dict)
        self.assertEqual(len(removed_bad_reads), 3)
        self.assertEqual(len(removed_bad_reads['contains_abnormal_reads.tsv']['AASeq']), 4) # We have not removed the long or short sequences.

    def test_remove_bad_len(self):
        ''' Tests that sequences that are below and above the biological threshold are removed. '''
        repertoires = self.cleaner.collect_files()
        removed_bad_reads = self.cleaner.remove_bad_reads(repertoires)
        removed_bad_len = self.cleaner.assert_size(removed_bad_reads)
        self.assertEqual(len(removed_bad_len), 3)
        self.assertEqual(len(removed_bad_len['contains_abnormal_reads.tsv']['AASeq']), 2) # We have removed the bad reads AND long and short sequences.
        self.assertEqual(len(removed_bad_len['123_no_errors.tsv']['AASeq']), 3)
        self.assertEqual(len(removed_bad_len['345_no_errors.tsv']['AASeq']), 3)

    def test_not_supported_disease(self):
        ''' Tests that an error is raised if the disease is not supported (meaning the model has not been trained on the data). '''
        repertoires = self.cleaner.collect_files()
        removed_bad_reads = self.cleaner.remove_bad_reads(repertoires)
        removed_bad_len = self.cleaner.assert_size(removed_bad_reads)
        removed_bad_len['This_is_not_a_disease'] = None
        self.assertRaises(DiseaseNotSupportedError, self.cleaner.join_files, removed_bad_len)


class TestExploratoryAnalysis(unittest.TestCase):
    ''' Tests the EDA Class. '''
    
    def setUp(self):
        handler = DataHandler(collection_path=Path('mock_data_for_tests/eda'),
                              save_path=Path('This_is_a_mock_path')) # Avoids the warning message.
        self.eda = ExploratoryAnalyzer(handler)
    
    def test_count_seqs_per_disease(self):
        ''' Tests the 'count_seqs_per_disease()' function. '''

        collected_data = self.eda.collect_repertoires()
        self.eda.count_sequences_per_disease(collected_data)
        self.assertEqual(self.eda.sequence_counts['test_1'], 3)
        self.assertEqual(self.eda.sequence_counts['test_2'], 5)
        self.assertEqual(self.eda.sequence_counts['test_3'], 2)

    def test_sequence_length_metrics(self):
        ''' Tests the function that counts the sequcence length metrics.'''
        collected_data = self.eda.collect_repertoires()
        self.eda.count_sequences_per_disease(collected_data)
        self.eda.count_sequence_length_metrics(collected_data)
        self.assertEqual(list(self.eda.seq_len_metrics['test_1']['len_counts']), [1.0])
        self.assertEqual(list(self.eda.seq_len_metrics['test_2']['len_counts']), [0.4, 0.2, 0.2, 0.2])
        self.assertEqual(list(self.eda.seq_len_metrics['test_3']['len_counts']), [0.5, 0.5])

    def test_gene_usage(self):
        ''' Tests the function that counts the gene usage. '''

        collected_data = self.eda.collect_repertoires()
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
        ''' Tests the function that computes the Jaccard Index. '''
        collected_data = self.eda.collect_repertoires()
        self.eda.compute_jaccard_index(collected_data)
        self.assertEqual(self.eda.jaccard_index['test_1'],  {'test_2': 0.0, 'test_3': 0.0})
        self.assertEqual(self.eda.jaccard_index['test_2'], {'test_1': 0.0, 'test_3': 0.16666666666666666})
        self.assertEqual(self.eda.jaccard_index['test_3'], {'test_1': 0.0, 'test_2': 0.16666666666666666})

    def test_aa_occurrences(self):
        ''' Tests the function that counts the amino acid occurrences.'''
        collected_data = self.eda.collect_repertoires()
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


# class TestKmerPreprocessing(unittest.TestCase):
#     ''' Tests the Preprocessing Class. '''

#     def setUp(self):
#         pass 



unittest.main()