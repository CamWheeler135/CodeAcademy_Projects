''' File that contains the code for cleaning data. '''

# Module Imports. 
from pathlib import Path, PurePath
import pandas as pd
import numpy as np
from tqdm import tqdm


# Local Imports
from base_logger import logging

class DataCleaner():
    '''
    This class is for cleaning and standardizing the data for use in the ML pipeline. 
    - Collects the raw data from the 'raw_data' file.
    - Cleans the data, removing sequences based on length and read quality.
    - Names files appropriately for clear analysis (file names can be set to remain the same if required).
    - Returns cleaned CSV files to 'cleaned_files'.
    - Outputs general dataset information to the terminal. 
    '''

    # Class variables.
    logger = logging.getLogger('Data Cleaner')
    number_of_data_files = 0
    supported_disease_types = ['colorectal', 'hematological_cancer', 'breast',
                    'lymphoblastic_leukemia', 'myeloid_leukemia', 'liver',
                    'lung', 'glioblastoma', 'gastric', 'esophageal']
    
    def __init__(self, input_path=Path('data/raw_files'), output_path=Path('data/cleaned_files')):
        self.input_path = input_path
        self.output_path = output_path
        self.min_seq_size = 10
        self.max_seq_size = 24
        self.special_characters = [r'\+', r'\*', r'_']
        self.column_order = ['AASeq', 'Vregion', 'Jregion', 'Dregion']

    def __str__(self):
        return "Object for general dataset cleaning, standard input --> data/raw_files, standard output --> data/cleaned_files. "

    @classmethod
    def collect_general_info(self):
        ''' Outputs general dataset information. '''
        if self.number_of_data_files == 0:
            self.logger.warning(f"There are currently no DataCleaner objects or the DataCleaner has not been called to clean any data.")
        else:
            self.logger.info("General Information about the dataset.")
            self.logger.info(f'The number of files that have been cleaned are {DataCleaner.number_of_data_files}')
            self.logger.info(f'The diseases that are currently accepted in the pipeline are {DataCleaner.supported_disease_types}')

    def collect_data(self) -> dict:
        ''' 
        Collect data from original files.
        Assumes that the file is from the TCRdb, support for other sources will be added later. 
        '''
        self.logger.info('Collecting files.')

        collection = {}
        for file in self.input_path.iterdir():
            self.logger.debug(f"Collecting {file.name}")

            # Alter the class variables. 
            DataCleaner.number_of_data_files += 1

            # Complete project files do not have certain metadata.
            if 'PRJ' in file.name:
                collection[file.name] = pd.read_table(file, header=0).drop(columns=['RunId', 'cloneFraction'])
            # Single files contain metadata, need to skip over. 
            else:
                collection[file.name] = pd.read_table(file, header=1)
            
            # Certain file column order are not the same, this code orders them correctly. 
            if not collection[file.name].columns.equals(self.column_order):
                self.logger.debug(f"Column order on {file.name} is not standard, changes have been made.")
                collection[file.name] = collection[file.name].reindex(columns=self.column_order)
            
        self.logger.info("All files collected.")
        
        return collection

    def assert_size(self, dataset:dict) -> dict:
        ''' 
        Remove sequences that exceed biological norms.
        - Short <= 10 
        - Long >= 24
        Reference: 
        Beshnova D, Ye J, Onabolu O, Moon B, Zheng W, Fu YX, Brugarolas J, Lea J, Li B. 
        De novo prediction of cancer-associated T cell receptors for noninvasive cancer detection. 
        '''
        self.logger.info("Searching for sequences that are outside of min/max range.")

        files = dataset.keys()
        for file in files:
            self.logger.debug(f"Handling {file}")
            self.logger.debug(f"The length of {file}'s data-frame before length removal = {len(dataset[file].index)}")

            # Remove sequences based on length parameters
            dataset[file] = dataset[file][(dataset[file]['AASeq'].str.len() >= self.min_seq_size) 
                                        & (dataset[file]['AASeq'].str.len() <= self.max_seq_size)]
            
            self.logger.debug(f"The length of {file}'s data-frame after length removal = {len(dataset[file].index)}")
        
        return dataset 

    def remove_bad_reads(self, dataset:dict, special_characters=None) -> dict:
        ''' 
        Remove sequences that contain non-alphabetical characters.
        All sequences from the TCRdb have been cleaned, but this will be needed for raw sequences. 
        Characters for removal:
        - '+'
        - '*'
        - '_'
        '''
        self.logger.info("Searching for sequences that contain bad reads.")

        files = dataset.keys()
        for file in files:
            self.logger.debug(f"Handling {file}.")
            for character in self.special_characters:
                self.logger.debug(f"Handling {character}")
                # Collect the index of the strings that contain errors. 
                idx = dataset[file].index[dataset[file]['AASeq'].str.contains(character)]
                if len(idx) == 0:
                    self.logger.debug("No special characters found.")
                    continue
                else:
                    # Drop the rows that contain the errors.
                    dataset[file].drop(idx, inplace=True)

        return dataset

    def join_files(self, cleaned_data) -> dict:
        '''
        Iterates through the dataset and aggregates the CDR3 sequences from the sample cancers into one data-frame.
        The list of cancers that are supported are held in the class variable 'supported_disease_types'. 
        '''
        self.logger.info("Creating final data frames.")

        # Take each cancer type and make one large data frame containing the data. 
        for cancer_type in DataCleaner.supported_disease_types:
            self.logger.debug(f"Handling {cancer_type} cancer files.")
            cancer_type_files = [file for file in cleaned_data.keys() if cancer_type in file]

            if len(cancer_type_files) <= 1:
                save_path = self.name_files(cancer_type)
                cleaned_data[cancer_type_files[0]].to_csv(save_path)
                self.logger.debug(f"Saved {cancer_type} data frame to {save_path}.")
            else:
            # Create the data frame for that cancer type and append all the data.)
                # print([cleaned_data[file].head() for file in cancer_type_files])
                cancer_type_df = pd.concat([cleaned_data[file] for file in cancer_type_files], ignore_index=True)
                save_path = self.name_files(cancer_type)
                cancer_type_df.to_csv(save_path)

    def name_files(self, cancer_type:str):
        ''' Gives files appropriate names according to disease. '''
        return PurePath('data/cleaned_files/' + cancer_type + '.csv')
        
    def complete_clean(self):
        ''' Performs a complete collection, cleaning and naming of dataset inside of input path. '''
        dataset = self.collect_data()
        bad_reads_cleaned = self.remove_bad_reads(dataset)
        length_cleaned = self.assert_size(bad_reads_cleaned)
        self.join_files(length_cleaned)
        self.collect_general_info()