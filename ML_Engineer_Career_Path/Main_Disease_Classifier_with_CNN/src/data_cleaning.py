''' File that contains the code for cleaning data. '''

# Module Imports. 
from pathlib import Path
import pandas as pd
import os


# Local Imports
from src.base_logger import logging
from src.errors import FileNotFoundError, DiseaseNotSupportedError

class DataHandler():
    '''
    This class is for handling files, primarily tasked with loading and saving.
    - Will load in the raw data from the TCRDB and read it into a dictionary of data frames.
    - Will load in cleaned data and read it into a dictionary of data frames. 
    - Will save data as CSV files.
    '''

    # Function that loads each of the CSV files for analysis.
    def __init__(self, collection_path:Path, save_path:Path=None):
        self.logger = logging.getLogger('Data Handler')
        self.collection_path = collection_path
        self.save_path = save_path
        if self.save_path == None:
            self.logger.warning('No save path given.')

    def collect_raw_tcrdb_files(self) -> dict:
        '''
        Collects the TCRDB data from the collection path file and returns a dictionary of the data.
        Assumes certain file properties:
        - The file is from the TCRdb.
        - Files are in the .tsv format.
        - File names contain the disease type. 
        Output {sample_name: pd.DataFrame}
        DataFrame columns (in order) = ['AASeq', 'Vregion', 'Dregion', 'Jregion']
        '''

        repertoire_dict = {}

        if len(os.listdir(self.collection_path)) == 0:
            raise FileNotFoundError

        for sample in self.collection_path.iterdir():
            self.logger.debug(f"Loading {sample.name} into dictionary.")

            # Handle TCRDB project files.
            if 'PRJ' in sample.name:
                repertoire_dict[sample.name] = pd.read_table(sample, header=0, 
                                                        usecols=['AASeq', 'Vregion', 'Dregion', 'Jregion', 'cloneFraction'])
            '''TODO'''
            # Support for other file types will be added later.

        '''TODO'''
        # We need to catch an error here, if the correct columns are not found, then we need to know how to process them. 

        self.logger.info("All files collected.")
        
        return repertoire_dict

    def collect_single_csv_file(self) -> pd.DataFrame:
        ''' Collects a single file and returns it as a data frame.'''
        return pd.read_csv(self.collection_path)

    def collect_csv_files(self, multi:bool=True) -> dict:
        '''Collects the data from the collection path file and returns a dictionary of the data.'''
        # Returns the file name (stripping the csv tag) with a data frame of the sequences.
        '''TODO'''
        # Bug currently does not remove teh CSV from the file name.
        return {os.path.splitext(file.name)[0]: pd.read_csv(file) for file in self.collection_path.iterdir()}

    def save_as_csv(self, data:pd.DataFrame):
        ''' Saves the data as a csv file. '''
        data.to_csv(self.save_path, index=False)
    
    def save_multi_as_csv(self, data:dict):
        ''' Iterates through a dictionary and saves the data as a csv. '''
        for file_name, data_df in data.items():
            data_df.to_csv(self.save_path / f'{file_name}.csv', index=False)

    '''TODO'''
    # Add support for error handling when collecting and saving files. 
    # Such as creating directories if we need to, handling empty files etc. 


class DataCleaner():
    '''
    This class is for cleaning and standardizing the data for use in the ML pipeline. 
    - Cleans the data, removing sequences based on length and read quality.
    - Names files appropriately for clear analysis (file names can be set to remain the same if required).
    - Outputs general dataset information to the terminal. 
    '''

    # Class variables.
    number_of_data_files = 0
    supported_disease_types = ['breast', 'lymphoblastic_leukemia',
                                'myeloid_leukemia', 'liver', 'lung', 'glioblastoma', 
                                'esophageal']
    
    def __init__(self, handler: DataHandler):
        self.logger = logging.getLogger('Data Cleaner')
        self.handler = handler
        self.min_seq_size = 10
        self.max_seq_size = 24
        self.special_characters = [r'\+', r'\*', r'_']

    def __str__(self):
        return f"Data Cleaner Object. \n\tMin Sequence Size: {self.min_seq_size} \n\tMax Sequence Size: {self.max_seq_size} \n\tSpecial Characters: {self.special_characters}"
    
    def collect_files(self) -> dict:
        ''' Collects the files from the data/raw_files directory using the DataHandler class. '''
        return self.handler.collect_raw_tcrdb_files()
    
    def remove_special_characters(self, dataset:dict, special_characters:list=None) -> dict:
        ''' 
        Remove sequences that contain non-alphabetical characters.
        All sequences from the TCRdb have been cleaned, but this function will be helpful with raw sequences being added to the project. 
        Characters for removal:
        - '+'
        - '*'
        - '_'
        '''
        self.logger.info("Searching for sequences that contain bad reads.")

        # Option to pass in specific special characters.
        if special_characters == None:
            special_characters = self.special_characters
        
        for file in dataset.keys():
            self.logger.debug(f"Handling {file}.")
            for character in special_characters:
                self.logger.debug(f"Handling {character}")
                # Collect the index of the strings that contain errors. 
                idx = dataset[file].index[dataset[file]['AASeq'].str.contains(character)]
                if len(idx) == 0:
                    self.logger.debug("No special characters found.")
                else:
                    # Drop the rows that contain the errors.
                    dataset[file] = dataset[file].drop(idx)

        return dataset

    def assert_size(self, dataset:dict) -> dict:
        ''' 
        Remove sequences that exceed biological norms.
        - Short <= 10 
        - Long >= 24
        Reference: 
        Beshnova D, et al.
        De novo prediction of cancer-associated T cell receptors for noninvasive cancer detection. 
        '''
        self.logger.info("Searching for sequences that are outside of min/max range.")

        files = dataset.keys()
        for file in files:
            self.logger.debug(f"The length of {file}'s data-frame before length removal = {len(dataset[file].index)}")

            # Strip the data of any possible whitespace.
            dataset[file]['AASeq'] = dataset[file]['AASeq'].str.strip()
            # Remove sequences based on length parameters
            dataset[file] = dataset[file][(dataset[file]['AASeq'].str.len() >= self.min_seq_size) 
                                        & (dataset[file]['AASeq'].str.len() <= self.max_seq_size)]
            self.logger.debug(f"The length of {file}'s data-frame after length removal = {len(dataset[file].index)}")
        
        return dataset 

    def remove_incomplete_seqs(self, dataset:dict) -> dict:
        '''
        Remove sequences that do not comply with IMGT standards.
        - Do not start with cysteine C. 
        - Do not end with phenylalanine F.
        Reference: 
        Xu, Ying et al. 
        DeepLION: Deep Multi-Instance Learning Improves the Prediction of Cancer-Associated T Cell Receptors 
        for Accurate Cancer Detection.
        '''
        self.logger.info("Searching for sequences that do not comply with IMGT standards.")

        for repertoire in dataset.keys():
            self.logger.debug(f"The length of the {repertoire} data-frame before IMGT removal = {len(dataset[repertoire].index)}")
            dataset[repertoire] = dataset[repertoire][(dataset[repertoire]['AASeq'].str.startswith('C'))
                                                      & (dataset[repertoire]['AASeq'].str.endswith('F'))]
            
            self.logger.debug(dataset[repertoire].head())
            self.logger.debug(f"The length of the {repertoire} data-frame after IMGT removal = {len(dataset[repertoire].index)}")

        return dataset

    def check_for_support(self, disease_to_check:str) -> bool:
        ''' 
        Checks if the disease is supported by the pipeline. 
        Returns the disease if it is supported.
        If not, raise an error. 
        '''
        for disease in DataCleaner.supported_disease_types:
            if disease in disease_to_check:
                return disease
        raise DiseaseNotSupportedError
    
    def name_files(self, cancer_type:str):
        ''' Gives files appropriate names according to disease. '''
        return Path('data/cleaned_files/' + cancer_type + '.csv')
    
    
    def join_files(self, cleaned_data):
        '''
        Iterates through the dataset and aggregates the CDR3 sequences from the sample cancers into one data-frame and saves them.
        The list of cancers that are supported are held in the class variable 'supported_disease_types'. 
        The complete data frames are saved in the handlers output_path.
        '''
        self.logger.info("Creating final data frames.")

        # Take each cancer type and make one large data frame containing the data.  
        for cancer_type in cleaned_data.keys():

            # Check if the cancer type is supported by the pipeline.
            cancer_to_concat = self.check_for_support(cancer_type)

            '''TODO'''
            # Add support to catch error. 

            self.logger.debug(f"Building {cancer_to_concat} data frame.")
            cancer_type_files = [file for file in cleaned_data.keys() if cancer_to_concat in file]
            self.logger.debug(f"Files to be joined: {cancer_type_files}")

            if len(cancer_type_files) == 1:
                save_path = self.name_files(cancer_to_concat)
                single_file = cancer_type_files[0]
                # Does not implement the handlers save function.
                cleaned_data[single_file].to_csv(save_path, index=False)
                self.logger.debug(f"Saved {cancer_type} data frame to {save_path}.")
            else:
                # Create the data frame for that cancer type and append all the data.
                cancer_type_df = pd.concat([cleaned_data[file] for file in cancer_type_files], ignore_index=True)
                save_path = self.name_files(cancer_to_concat)
                # Does not implement the handlers save function.
                cancer_type_df.to_csv(save_path, index=False)
                self.logger.debug(f"Saved {cancer_type} data frame to {save_path}.")
        
        self.logger.info("Final data frames saved.")
        
    def complete_clean(self):
        ''' Performs a complete collection, cleaning and naming of dataset inside of input path. '''
        dataset = self.collect_files()
        bad_reads_cleaned = self.remove_special_characters(dataset)
        length_cleaned = self.assert_size(bad_reads_cleaned)
        complete_seqs = self.remove_incomplete_seqs(length_cleaned)
        self.join_files(complete_seqs)
        self.logger.info("Data cleaning complete.")