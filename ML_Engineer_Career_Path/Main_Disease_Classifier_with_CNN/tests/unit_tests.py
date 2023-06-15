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


class TestKmerPreprocessing(unittest.TestCase):
    ''' Tests the Preprocessing Class. '''

    def setUp(self):
        kmer_handler = DataHandler(collection_path=Path('mock_data_for_tests/preprocessing/kmer_data'),
                                   save_path=Path('This_is_a_mock_path')) # Avoids the warning log. 
        self.kmer_preprocessor = KmerPreprocessor(kmer_handler, subsample_size=2, num_of_samples=2)

    def test_subampling(self):
        ''' Tests that the code will subsample the correct amount of sequences from each repertoire.'''
        repertoires = self.kmer_preprocessor.collect_repertoires()
        subsamples = self.kmer_preprocessor.sub_sample_repertoire(repertoires)
        self.assertEqual(len(subsamples), 2) # There should only be 2 datasets.
        self.assertEqual(len(subsamples['test_1']), 2) # There should be 2 sublists in value list.
        self.assertEqual(len(subsamples['test_1'][0]), 2) # There should be 2 sequences in each sublist.


    def test_target_creation(self):
        ''' Tests that the targets are created appropriately. '''
        repertoires = self.kmer_preprocessor.collect_repertoires()
        subsamples = self.kmer_preprocessor.sub_sample_repertoire(repertoires)
        target_vals = self.kmer_preprocessor.create_targets(subsamples)
        self.assertEqual(len(target_vals), 2) # There should only be 2 targets.
        self.assertEqual(target_vals['test_1'], 0) # The first target should be 0.
        self.assertEqual(target_vals['test_2'], 1) # The second target should be 1.


    # We need to test that we are creating the kmer counts correctly.
    def test_kmer_counts(self):
        ''' Tests that the correct number of kmers are counted for a sequence. '''
        # Hard code seqs so we have a known output.
        sequences = {'test_1': [['CASSLGGDFDTEAFF']], 
                     'test_2': [['CASSLGGDFDTEAFF']]}
        
        # Hard code the kmer dict to test the output, if this changes, we are not counting the kmers correctly.
        seq1_kmer_count = {'AA': 0, 'AR': 0, 'AN': 0, 'AD': 0, 'AC': 0, 'AE': 0, 'AQ': 0, 'AG': 0, 'AH': 0, 'AI': 0, 'AL': 0, 'AK': 0, 
                            'AM': 0, 'AF': 1, 'AP': 0, 'AS': 1, 'AT': 0, 'AW': 0, 'AY': 0, 'AV': 0, 'RA': 0, 'RR': 0, 'RN': 0, 'RD': 0, 
                            'RC': 0, 'RE': 0, 'RQ': 0, 'RG': 0, 'RH': 0, 'RI': 0, 'RL': 0, 'RK': 0, 'RM': 0, 'RF': 0, 'RP': 0, 'RS': 0, 
                            'RT': 0, 'RW': 0, 'RY': 0, 'RV': 0, 'NA': 0, 'NR': 0, 'NN': 0, 'ND': 0, 'NC': 0, 'NE': 0, 'NQ': 0, 'NG': 0, 
                            'NH': 0, 'NI': 0, 'NL': 0, 'NK': 0, 'NM': 0, 'NF': 0, 'NP': 0, 'NS': 0, 'NT': 0, 'NW': 0, 'NY': 0, 'NV': 0, 
                            'DA': 0, 'DR': 0, 'DN': 0, 'DD': 0, 'DC': 0, 'DE': 0, 'DQ': 0, 'DG': 0, 'DH': 0, 'DI': 0, 'DL': 0, 'DK': 0, 
                            'DM': 0, 'DF': 1, 'DP': 0, 'DS': 0, 'DT': 1, 'DW': 0, 'DY': 0, 'DV': 0, 'CA': 1, 'CR': 0, 'CN': 0, 'CD': 0, 
                            'CC': 0, 'CE': 0, 'CQ': 0, 'CG': 0, 'CH': 0, 'CI': 0, 'CL': 0, 'CK': 0, 'CM': 0, 'CF': 0, 'CP': 0, 'CS': 0, 
                            'CT': 0, 'CW': 0, 'CY': 0, 'CV': 0, 'EA': 1, 'ER': 0, 'EN': 0, 'ED': 0, 'EC': 0, 'EE': 0, 'EQ': 0, 'EG': 0, 
                            'EH': 0, 'EI': 0, 'EL': 0, 'EK': 0, 'EM': 0, 'EF': 0, 'EP': 0, 'ES': 0, 'ET': 0, 'EW': 0, 'EY': 0, 'EV': 0, 
                            'QA': 0, 'QR': 0, 'QN': 0, 'QD': 0, 'QC': 0, 'QE': 0, 'QQ': 0, 'QG': 0, 'QH': 0, 'QI': 0, 'QL': 0, 'QK': 0,
                            'QM': 0, 'QF': 0, 'QP': 0, 'QS': 0, 'QT': 0, 'QW': 0, 'QY': 0, 'QV': 0, 'GA': 0, 'GR': 0, 'GN': 0, 'GD': 1, 
                            'GC': 0, 'GE': 0, 'GQ': 0, 'GG': 1, 'GH': 0, 'GI': 0, 'GL': 0, 'GK': 0, 'GM': 0, 'GF': 0, 'GP': 0, 'GS': 0, 
                            'GT': 0, 'GW': 0, 'GY': 0, 'GV': 0, 'HA': 0, 'HR': 0, 'HN': 0, 'HD': 0, 'HC': 0, 'HE': 0, 'HQ': 0, 'HG': 0, 
                            'HH': 0, 'HI': 0, 'HL': 0, 'HK': 0, 'HM': 0, 'HF': 0, 'HP': 0, 'HS': 0, 'HT': 0, 'HW': 0, 'HY': 0, 'HV': 0, 
                            'IA': 0, 'IR': 0, 'IN': 0, 'ID': 0, 'IC': 0, 'IE': 0, 'IQ': 0, 'IG': 0, 'IH': 0, 'II': 0, 'IL': 0, 'IK': 0, 
                            'IM': 0, 'IF': 0, 'IP': 0, 'IS': 0, 'IT': 0, 'IW': 0, 'IY': 0, 'IV': 0, 'LA': 0, 'LR': 0, 'LN': 0, 'LD': 0, 
                            'LC': 0, 'LE': 0, 'LQ': 0, 'LG': 1, 'LH': 0, 'LI': 0, 'LL': 0, 'LK': 0, 'LM': 0, 'LF': 0, 'LP': 0, 'LS': 0, 
                            'LT': 0, 'LW': 0, 'LY': 0, 'LV': 0, 'KA': 0, 'KR': 0, 'KN': 0, 'KD': 0, 'KC': 0, 'KE': 0, 'KQ': 0, 'KG': 0, 
                            'KH': 0, 'KI': 0, 'KL': 0, 'KK': 0, 'KM': 0, 'KF': 0, 'KP': 0, 'KS': 0, 'KT': 0, 'KW': 0, 'KY': 0, 'KV': 0, 
                            'MA': 0, 'MR': 0, 'MN': 0, 'MD': 0, 'MC': 0, 'ME': 0, 'MQ': 0, 'MG': 0, 'MH': 0, 'MI': 0, 'ML': 0, 'MK': 0, 
                            'MM': 0, 'MF': 0, 'MP': 0, 'MS': 0, 'MT': 0, 'MW': 0, 'MY': 0, 'MV': 0, 'FA': 0, 'FR': 0, 'FN': 0, 'FD': 1, 
                            'FC': 0, 'FE': 0, 'FQ': 0, 'FG': 0, 'FH': 0, 'FI': 0, 'FL': 0, 'FK': 0, 'FM': 0, 'FF': 1, 'FP': 0, 'FS': 0, 
                            'FT': 0, 'FW': 0, 'FY': 0, 'FV': 0, 'PA': 0, 'PR': 0, 'PN': 0, 'PD': 0, 'PC': 0, 'PE': 0, 'PQ': 0, 'PG': 0, 
                            'PH': 0, 'PI': 0, 'PL': 0, 'PK': 0, 'PM': 0, 'PF': 0, 'PP': 0, 'PS': 0, 'PT': 0, 'PW': 0, 'PY': 0, 'PV': 0, 
                            'SA': 0, 'SR': 0, 'SN': 0, 'SD': 0, 'SC': 0, 'SE': 0, 'SQ': 0, 'SG': 0, 'SH': 0, 'SI': 0, 'SL': 1, 'SK': 0, 
                            'SM': 0, 'SF': 0, 'SP': 0, 'SS': 1, 'ST': 0, 'SW': 0, 'SY': 0, 'SV': 0, 'TA': 0, 'TR': 0, 'TN': 0, 'TD': 0, 
                            'TC': 0, 'TE': 1, 'TQ': 0, 'TG': 0, 'TH': 0, 'TI': 0, 'TL': 0, 'TK': 0, 'TM': 0, 'TF': 0, 'TP': 0, 'TS': 0, 
                            'TT': 0, 'TW': 0, 'TY': 0, 'TV': 0, 'WA': 0, 'WR': 0, 'WN': 0, 'WD': 0, 'WC': 0, 'WE': 0, 'WQ': 0, 'WG': 0, 
                            'WH': 0, 'WI': 0, 'WL': 0, 'WK': 0, 'WM': 0, 'WF': 0, 'WP': 0, 'WS': 0, 'WT': 0, 'WW': 0, 'WY': 0, 'WV': 0, 
                            'YA': 0, 'YR': 0, 'YN': 0, 'YD': 0, 'YC': 0, 'YE': 0, 'YQ': 0, 'YG': 0, 'YH': 0, 'YI': 0, 'YL': 0, 'YK': 0, 
                            'YM': 0, 'YF': 0, 'YP': 0, 'YS': 0, 'YT': 0, 'YW': 0, 'YY': 0, 'YV': 0, 'VA': 0, 'VR': 0, 'VN': 0, 'VD': 0, 
                            'VC': 0, 'VE': 0, 'VQ': 0, 'VG': 0, 'VH': 0, 'VI': 0, 'VL': 0, 'VK': 0, 'VM': 0, 'VF': 0, 'VP': 0, 'VS': 0, 
                            'VT': 0, 'VW': 0, 'VY': 0, 'VV': 0, 'target': 0}

        seq2_kmer_count = {'AA': 0, 'AR': 0, 'AN': 0, 'AD': 0, 'AC': 0, 'AE': 0, 'AQ': 0, 'AG': 0, 'AH': 0, 'AI': 0, 'AL': 0, 'AK': 0, 
                           'AM': 0, 'AF': 1, 'AP': 0, 'AS': 1, 'AT': 0, 'AW': 0, 'AY': 0, 'AV': 0, 'RA': 0, 'RR': 0, 'RN': 0, 'RD': 0, 
                           'RC': 0, 'RE': 0, 'RQ': 0, 'RG': 0, 'RH': 0, 'RI': 0, 'RL': 0, 'RK': 0, 'RM': 0, 'RF': 0, 'RP': 0, 'RS': 0, 
                           'RT': 0, 'RW': 0, 'RY': 0, 'RV': 0, 'NA': 0, 'NR': 0, 'NN': 0, 'ND': 0, 'NC': 0, 'NE': 0, 'NQ': 0, 'NG': 0, 
                           'NH': 0, 'NI': 0, 'NL': 0, 'NK': 0, 'NM': 0, 'NF': 0, 'NP': 0, 'NS': 0, 'NT': 0, 'NW': 0, 'NY': 0, 'NV': 0, 
                           'DA': 0, 'DR': 0, 'DN': 0, 'DD': 0, 'DC': 0, 'DE': 0, 'DQ': 0, 'DG': 0, 'DH': 0, 'DI': 0, 'DL': 0, 'DK': 0, 
                           'DM': 0, 'DF': 1, 'DP': 0, 'DS': 0, 'DT': 1, 'DW': 0, 'DY': 0, 'DV': 0, 'CA': 1, 'CR': 0, 'CN': 0, 'CD': 0, 
                           'CC': 0, 'CE': 0, 'CQ': 0, 'CG': 0, 'CH': 0, 'CI': 0, 'CL': 0, 'CK': 0, 'CM': 0, 'CF': 0, 'CP': 0, 'CS': 0, 
                           'CT': 0, 'CW': 0, 'CY': 0, 'CV': 0, 'EA': 1, 'ER': 0, 'EN': 0, 'ED': 0, 'EC': 0, 'EE': 0, 'EQ': 0, 'EG': 0, 
                           'EH': 0, 'EI': 0, 'EL': 0, 'EK': 0, 'EM': 0, 'EF': 0, 'EP': 0, 'ES': 0, 'ET': 0, 'EW': 0, 'EY': 0, 'EV': 0, 
                           'QA': 0, 'QR': 0, 'QN': 0, 'QD': 0, 'QC': 0, 'QE': 0, 'QQ': 0, 'QG': 0, 'QH': 0, 'QI': 0, 'QL': 0, 'QK': 0, 
                           'QM': 0, 'QF': 0, 'QP': 0, 'QS': 0, 'QT': 0, 'QW': 0, 'QY': 0, 'QV': 0, 'GA': 0, 'GR': 0, 'GN': 0, 'GD': 1, 
                           'GC': 0, 'GE': 0, 'GQ': 0, 'GG': 1, 'GH': 0, 'GI': 0, 'GL': 0, 'GK': 0, 'GM': 0, 'GF': 0, 'GP': 0, 'GS': 0, 
                           'GT': 0, 'GW': 0, 'GY': 0, 'GV': 0, 'HA': 0, 'HR': 0, 'HN': 0, 'HD': 0, 'HC': 0, 'HE': 0, 'HQ': 0, 'HG': 0, 
                           'HH': 0, 'HI': 0, 'HL': 0, 'HK': 0, 'HM': 0, 'HF': 0, 'HP': 0, 'HS': 0, 'HT': 0, 'HW': 0, 'HY': 0, 'HV': 0, 
                           'IA': 0, 'IR': 0, 'IN': 0, 'ID': 0, 'IC': 0, 'IE': 0, 'IQ': 0, 'IG': 0, 'IH': 0, 'II': 0, 'IL': 0, 'IK': 0, 
                           'IM': 0, 'IF': 0, 'IP': 0, 'IS': 0, 'IT': 0, 'IW': 0, 'IY': 0, 'IV': 0, 'LA': 0, 'LR': 0, 'LN': 0, 'LD': 0, 
                           'LC': 0, 'LE': 0, 'LQ': 0, 'LG': 1, 'LH': 0, 'LI': 0, 'LL': 0, 'LK': 0, 'LM': 0, 'LF': 0, 'LP': 0, 'LS': 0, 
                           'LT': 0, 'LW': 0, 'LY': 0, 'LV': 0, 'KA': 0, 'KR': 0, 'KN': 0, 'KD': 0, 'KC': 0, 'KE': 0, 'KQ': 0, 'KG': 0, 
                           'KH': 0, 'KI': 0, 'KL': 0, 'KK': 0, 'KM': 0, 'KF': 0, 'KP': 0, 'KS': 0, 'KT': 0, 'KW': 0, 'KY': 0, 'KV': 0, 
                           'MA': 0, 'MR': 0, 'MN': 0, 'MD': 0, 'MC': 0, 'ME': 0, 'MQ': 0, 'MG': 0, 'MH': 0, 'MI': 0, 'ML': 0, 'MK': 0, 
                           'MM': 0, 'MF': 0, 'MP': 0, 'MS': 0, 'MT': 0, 'MW': 0, 'MY': 0, 'MV': 0, 'FA': 0, 'FR': 0, 'FN': 0, 'FD': 1, 
                           'FC': 0, 'FE': 0, 'FQ': 0, 'FG': 0, 'FH': 0, 'FI': 0, 'FL': 0, 'FK': 0, 'FM': 0, 'FF': 1, 'FP': 0, 'FS': 0, 
                           'FT': 0, 'FW': 0, 'FY': 0, 'FV': 0, 'PA': 0, 'PR': 0, 'PN': 0, 'PD': 0, 'PC': 0, 'PE': 0, 'PQ': 0, 'PG': 0, 
                           'PH': 0, 'PI': 0, 'PL': 0, 'PK': 0, 'PM': 0, 'PF': 0, 'PP': 0, 'PS': 0, 'PT': 0, 'PW': 0, 'PY': 0, 'PV': 0, 
                           'SA': 0, 'SR': 0, 'SN': 0, 'SD': 0, 'SC': 0, 'SE': 0, 'SQ': 0, 'SG': 0, 'SH': 0, 'SI': 0, 'SL': 1, 'SK': 0, 
                           'SM': 0, 'SF': 0, 'SP': 0, 'SS': 1, 'ST': 0, 'SW': 0, 'SY': 0, 'SV': 0, 'TA': 0, 'TR': 0, 'TN': 0, 'TD': 0, 
                           'TC': 0, 'TE': 1, 'TQ': 0, 'TG': 0, 'TH': 0, 'TI': 0, 'TL': 0, 'TK': 0, 'TM': 0, 'TF': 0, 'TP': 0, 'TS': 0, 
                           'TT': 0, 'TW': 0, 'TY': 0, 'TV': 0, 'WA': 0, 'WR': 0, 'WN': 0, 'WD': 0, 'WC': 0, 'WE': 0, 'WQ': 0, 'WG': 0, 
                           'WH': 0, 'WI': 0, 'WL': 0, 'WK': 0, 'WM': 0, 'WF': 0, 'WP': 0, 'WS': 0, 'WT': 0, 'WW': 0, 'WY': 0, 'WV': 0, 
                           'YA': 0, 'YR': 0, 'YN': 0, 'YD': 0, 'YC': 0, 'YE': 0, 'YQ': 0, 'YG': 0, 'YH': 0, 'YI': 0, 'YL': 0, 'YK': 0, 
                           'YM': 0, 'YF': 0, 'YP': 0, 'YS': 0, 'YT': 0, 'YW': 0, 'YY': 0, 'YV': 0, 'VA': 0, 'VR': 0, 'VN': 0, 'VD': 0, 
                           'VC': 0, 'VE': 0, 'VQ': 0, 'VG': 0, 'VH': 0, 'VI': 0, 'VL': 0, 'VK': 0, 'VM': 0, 'VF': 0, 'VP': 0, 'VS': 0, 
                           'VT': 0, 'VW': 0, 'VY': 0, 'VV': 0, 'target': 1}
        
        kmer_counts = self.kmer_preprocessor.collect_subsample_kmers(sequences)
        self.assertEqual(kmer_counts[0], seq1_kmer_count)
        self.assertEqual(kmer_counts[1], seq2_kmer_count)


unittest.main()