''' File that contains the code for exploratory data analysis and machine learning pre-processing. '''

# Module Imports.
from abc import ABC, abstractmethod
from pathlib import Path
import pandas as pd
import numpy as np

# Local Imports.
from src.base_logger import logging
from src.data_cleaning import DataLoader

class BasePreprocessor(ABC):

    @abstractmethod
    def sub_sample_repertoire(self):
        pass

    @abstractmethod
    def preprocess_data(self):
        pass

    @abstractmethod
    def save_data(self):
        pass


class UnsupervisedPreprocessor(BasePreprocessor, DataLoader):
    '''
    This class performs the preprocessing for the UnsupervisedVisualizer class.
    The preprocessing steps include:
    - Sub-sampling each repertoire into size n.
    - Counting the kmers present in each subsample. 
    - Saving the data frame as a CSV file in the data/processed/ directory.
    '''

    def __init__(self, subsample_size:int, kmer_size:int, collection_path=Path('data/cleaned_files')):
        self.logger = logging.getLogger('Unsupervised Preprocessor')
        super().__init__(collection_path) 
        self.subsample_size = subsample_size
        self.kmer_size = kmer_size
        self.kmers = {'''TODO'''}


    # Abstract method.
    def sub_sample_repertoire(self, data:dict) -> dict:
        ''' Take each repertoire and subsample it self.subsample_size times.'''
        pass

    def count_kmers(self, data:dict) -> dict:
        ''' For each subsample, count the number of kmers present.'''
        pass

    # Abstract method.
    def save_data(self):
        ''' Save the data as CSV file in the data/processed/ directory.'''
        pass

    # Abstract method.
    def preprocess_data(self):
        ''' 
        Preprocess the data into a format that can be used by the UnsupervisedVisualizer class.
        - Create the desired Kmers.
        - Take the sub-samples and count the kmers present in each subsample.
        '''
        pass


class CNNPreprocessor(BasePreprocessor, DataLoader):

    def __init__(self, subsample_size:int, collection_path=Path('data/cleaned_files')):
        self.logger = logging.getLogger('CNN Preprocessor')
        super().__init__(collection_path)
        pass

    # Abstract method.
    def sub_sample_repertoire(self):
        pass

    # Abstract method.
    def preprocess_data(self):
        pass

    # Abstract method.
    def save_data(self):
        pass