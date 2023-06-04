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
    Class will load and save the data into a dictionary for use in analysis.
    '''

    # Function that loads each of the CSV files for analysis.
    def __init__(self, collection_path: Path, save_path: Path):
        self.collection_path = collection_path
        self.save_path = save_path

    def collect_raw_files(self) -> dict:
        '''Collects the data from the collection path file and returns a dictionary of the data.'''
        return {os.path.splitext(file.name)[0]: pd.read_table(file) for file in self.collection_path.iterdir()}

    def collect_cleaned_files(self) -> dict:
        '''Collects the data from the collection path file and returns a dictionary of the data.'''
        # Returns the file name (stripping the csv tag) with a data frame of the sequences.
        return {os.path.splitext(file.name)[0]: pd.read_csv(file) for file in self.collection_path.iterdir()}
    
    def save_as_csv(self, data:pd.DataFrame):
        ''' Saves the data as a csv file. '''
        data.to_csv(self.save_path, index=False)

    '''TODO'''
    # Add support for error handling when collecting and saving files. 
    # Such as creating directories if we need to, handling empty files etc. 
    # This should help us take away the file handling from our classes that access and save data.



class DataCleaner():
    '''
    This class is for cleaning and standardizing the data for use in the ML pipeline. 
    - Cleans the data, removing sequences based on length and read quality.
    - Names files appropriately for clear analysis (file names can be set to remain the same if required).
    - Outputs general dataset information to the terminal. 
    '''

    # Class variables.
    number_of_data_files = 0
    supported_disease_types = ['colorectal', 'hematological_cancer', 'breast',
                    'lymphoblastic_leukemia', 'myeloid_leukemia', 'liver',
                    'lung', 'glioblastoma', 'esophageal']
    
    def __init__(self, handler: DataHandler):
        self.logger = logging.getLogger('Data Cleaner')
        self.handler = handler
        self.min_seq_size = 10
        self.max_seq_size = 24
        self.special_characters = [r'\+', r'\*', r'_']
        self.column_order = ['AASeq', 'Vregion', 'Jregion', 'Dregion']

    def __str__(self):
        return "Object for general dataset cleaning, standard input --> data/raw_files, standard output --> data/cleaned_files. "
    
    def collect_files(self) -> dict:
        return self.handler.collect_raw_files()

    def collect_data_for_analysis(self, repertoires: dict) -> dict:
        ''' 
        Remove un-needed columns, will reorder columns so each data frame is homogenous.
        There are several different files available on the TCRdb, this function will handle both single samples, and whole projects.
        The two types of files have different structures so need to be handled differently.
        Assumptions:
        - The file is from the TCRdb, support for other sources will be added later. 
        - Project files are labelled with the prefix PRJ.
        - File names contain the disease type. 
        '''
        self.logger.info('Collecting files.')

        if len(os.listdir(self.input_path)) == 0:
            raise FileNotFoundError

        collection = {}

        for sample, data in repertoires.items():

            pass

            '''TODO: Refactor this code to be more efficient.'''

        #     # Complete project files do not have certain metadata.
        #     if 'PRJ' in sample:
        #         collection[sample] = pd.read_table(file, header=0)
        #     # Single files contain metadata, need to skip over. 
        #     else:
        #         collection[sample] = pd.read_table(file, header=1)
             
        #     # Ensure that the data frames have the columns needed for analysis with EDA and ML. 
        #     if not set(self.column_order).issubset(collection[file.name].columns):
        #         self.logger.warning(f"{file.name} does not contain the correct columns for analysis, ensure columns are correct.")
        #         self.logger.warning(f"Removing {file.name} from analysis.")
        #         del collection[file.name]
        #         # Move onto the next file. 
        #         continue

        #     # Certain file column order are not the same, this code orders them correctly.
        #     if not collection[file.name].columns.equals(self.column_order):
        #         self.logger.debug(f"Column order on {file.name} is not standard, changes have been made.")
        #         collection[file.name] = collection[file.name].reindex(columns=self.column_order) # This will also drop any columns we don't need. 
            
        # self.logger.info("All files collected.")
        
        # return collection
    
    def remove_bad_reads(self, dataset:dict, special_characters:list=None) -> dict:
        ''' 
        Remove sequences that contain non-alphabetical characters.
        All sequences from the TCRdb have been cleaned, but this function will be helpful with raw sequences being added to the project. 
        Characters for removal:
        - '+'
        - '*'
        - '_'
        '''
        self.logger.info("Searching for sequences that contain bad reads.")

        files = dataset.keys()
        # Option to pass in specific special characters.
        if special_characters == None:
            special_characters = self.special_characters
        
        for file in files:
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
        Beshnova D, Ye J, Onabolu O, Moon B, Zheng W, Fu YX, Brugarolas J, Lea J, Li B. 
        De novo prediction of cancer-associated T cell receptors for noninvasive cancer detection. 
        '''
        self.logger.info("Searching for sequences that are outside of min/max range.")

        files = dataset.keys()
        for file in files:
            self.logger.debug(f"Handling {file}")
            self.logger.debug(f"The length of {file}'s data-frame before length removal = {len(dataset[file].index)}")

            # Strip the data of any possible whitespace.
            dataset[file]['AASeq'] = dataset[file]['AASeq'].str.strip()
            # Remove sequences based on length parameters
            dataset[file] = dataset[file][(dataset[file]['AASeq'].str.len() >= self.min_seq_size) 
                                        & (dataset[file]['AASeq'].str.len() <= self.max_seq_size)]
            
            self.logger.debug(f"The length of {file}'s data-frame after length removal = {len(dataset[file].index)}")
        
        return dataset 
    
    def join_files(self, cleaned_data):
        '''
        Iterates through the dataset and aggregates the CDR3 sequences from the sample cancers into one data-frame and saves them.
        The list of cancers that are supported are held in the class variable 'supported_disease_types'. 
        The complete data frames are saved in the instances output_path.
        '''
        self.logger.info("Creating final data frames.")

        # Take each cancer type and make one large data frame containing the data. 
        for cancer_type in DataCleaner.supported_disease_types:
            self.logger.info(f"Building {cancer_type} data frame.")
            self.logger.debug(f"Handling {cancer_type} cancer files.")
            cancer_type_files = [file for file in cleaned_data.keys() if cancer_type in file]

            if len(cancer_type_files) == 0:
                raise DiseaseNotSupportedError
            elif len(cancer_type_files) == 1:
                save_path = self.name_files(cancer_type)
                cleaned_data[cancer_type_files[0]].to_csv(save_path, index=False)
                self.logger.debug(f"Saved {cancer_type} data frame to {save_path}.")
            else:
            # Create the data frame for that cancer type and append all the data.)
                cancer_type_df = pd.concat([cleaned_data[file] for file in cancer_type_files], ignore_index=True)
                save_path = self.name_files(cancer_type)
                cancer_type_df.to_csv(save_path, index=False)
                self.logger.debug(f"Saved {cancer_type} data frame to {save_path}.")

    def name_files(self, cancer_type:str):
        ''' Gives files appropriate names according to disease. '''
        return Path('data/cleaned_files/' + cancer_type + '.csv')
        
    def complete_clean(self):
        ''' Performs a complete collection, cleaning and naming of dataset inside of input path. '''
        dataset = self.collect_data()
        bad_reads_cleaned = self.remove_bad_reads(dataset)
        length_cleaned = self.assert_size(bad_reads_cleaned)
        self.join_files(length_cleaned)
        self.collect_general_info()
        self.logger.info("Data cleaning complete.")