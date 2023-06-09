''' File that contains the code for exploratory data analysis and machine learning pre-processing. '''

# Module Imports.
from abc import ABC, abstractmethod
from pathlib import Path
import pandas as pd
import numpy as np
from itertools import product
import nltk

# Local Imports.
from src.base_logger import logging
from src.data_cleaning import DataHandler

# Skeletal must have methods for each preprocessing class.
class BasePreprocessor(ABC):

    @abstractmethod
    def collect_repertoires(self):
        pass

    @abstractmethod
    def sub_sample_repertoire(self):
        pass

    @abstractmethod
    def create_targets(self):
        pass

    @abstractmethod
    def save_data(self):
        pass

    @abstractmethod
    def complete_preprocess(self):
        pass


class KmerPreprocessor(BasePreprocessor):
    '''
    This class performs the preprocessing for the UnsupervisedVisualizer class.
    The preprocessing steps include:
    - Sub-sampling each repertoire x times into size n sub samples.
    - Counting the kmers of size m (default value m=2) present in each subsample. 
    - Constructing a data frame of kmers present in each subsample.
    - Receives a DataHandler object upon initialization to load save the data as a csv file.
    '''

    def __init__(
            self, handler: DataHandler, subsample_size:int=100, 
            num_of_samples:int=1000, kmer_size:int=2):
        
        self.logger = logging.getLogger('Unsupervised Preprocessor')
        self.subsample_size = subsample_size
        self.num_of_samples = num_of_samples
        self.kmer_size = kmer_size
        self.base_kmer_counter = {''.join(kmer): 0 for kmer in product('ARNDCEQGHILKMFPSTWYV', repeat=self.kmer_size)}
        self.handler = handler

    # Abstract method.
    def collect_repertoires(self) -> dict:
        return self.handler.collect_csv_files()

    # Abstract method.
    def sub_sample_repertoire(self, data:dict) -> dict:
        ''' 
        Take each repertoire and subsample it self.subsample_size times.
        Returns a dictionary of lists where each key is the disease name and each value is a list of lists.
        Each sub-list contains the sequences sampled from the repertoire.
        '''

        self.logger.info(f"Sub-sampling each repertoire into {self.num_of_samples} samples of size {self.subsample_size} seqs. ")
        repertoire_subsamples = {}
        for repertoire in data.keys():
            # Dictionary structure: {repertoire: [[subsample1], [subsample2], ...]}
            repertoire_subsamples[repertoire] = [data[repertoire]['AASeq'].sample(self.subsample_size, replace=True).to_list()
                                                 for i in range(self.num_of_samples)]
            self.logger.debug(f"The number of subsamples per reperotire = {len(repertoire_subsamples[repertoire])}")
            self.logger.debug(f"The number of sequences in each subsample = {len(repertoire_subsamples[repertoire][0])}")
        return repertoire_subsamples

    # Abstract method.
    def create_targets(self, data:dict) -> dict:
        ''' Creates a dictionary of the target values according to the dataset it has been given.'''
        return {disease: value for value, disease in enumerate(sorted(data.keys()))} 
    
    def count_kmers(self, sequences:list) -> dict:
        '''
        Iterates through each sequence in the sub-sample, counting the number of kmers present.
        Returns a dictionary of the kmers and their counts.
        '''
        # Construct the kmer dictionary.
        kmer_vals = self.base_kmer_counter.copy()
        # Iterate through the sequences in the subsample and increment the dictionary values.
        for sequence in sequences:
            kmer_iterator = nltk.ngrams(sequence, self.kmer_size)
            for kmer in kmer_iterator:
                kmer_vals[''.join(kmer)] += 1 
        return kmer_vals

    def collect_subsample_kmers(self, repertoire_subsamples:dict) -> list:
        ''' 
        For each subsample:
        - Count the number of kmers present.
        - Add the target value to each repertoire.
        - Append the results to a list.
        - Return the list.
        '''

        # Create the targets.
        target_vals = self.create_targets(repertoire_subsamples)
        self.logger.debug(f"Target values: {target_vals}")

        # Count the kmers and add the target to each dict.
        kmer_list = []
        for repertoire in repertoire_subsamples.keys():
            self.logger.info(f"Counting the kmers of size {self.kmer_size} present in each {repertoire} sub-sample.'")
            for subsample in repertoire_subsamples[repertoire]:
                subsample_kmers = self.count_kmers(subsample)
                # Adds the target value to the dictionary of kmers.
                subsample_kmers['target'] = target_vals[repertoire]
                kmer_list.append(subsample_kmers)

        return kmer_list

    # Abstract method.
    def construct_kmer_df(self, kmer_list:list) -> pd.DataFrame:
        ''' 
        Takes the list of kmer dictionaries.
        Forms a large data frame out of all the repertoires.
        Returns a data frame of kmers.
        '''
        kmer_df = pd.DataFrame(kmer_list)
        self.logger.debug(f"The shape of kmer data frame is {kmer_df.shape}")
        return kmer_df
    
    def save_data(self, kmer_df: pd.DataFrame):
        ''' Save the data frame as a CSV file in the handlers 'save_path' directory.'''
        self.handler.save_as_csv(kmer_df)

    # Abstract method.
    def complete_preprocess(self):
        ''' 
        Preprocess the data into a format that can be used by the UnsupervisedVisualizer class.
        - Create the desired Kmers.
        - Take the sub-samples and count the kmers present in each subsample.
        '''
        repertoires = self.collect_repertoires()
        subsamples = self.sub_sample_repertoire(repertoires)
        kmers = self.collect_subsample_kmers(subsamples)
        kmer_df = self.construct_kmer_df(kmers)
        self.save_data(kmer_df)


class CNNPreprocessor(BasePreprocessor):

    def __init__(
            self, handler: DataHandler, seqs_to_save:int=50000, 
            subsample_size:int=100, num_of_samples:int=10000):
        self.logger = logging.getLogger('CNN Preprocessor')
        self.handler = handler
        self.seqs_to_save = seqs_to_save
        self.subsample_size = subsample_size
        self.num_of_samples = num_of_samples
        self.encoding_matrix = None # Add the Beshova matrix here. 

    # Abstract method.
    def collect_repertoires(self):
        return self.handler.collect_csv_files()

    def save_top_n_seqs(self, data:dict) -> dict:
        '''
        Takes each repertoire, sorts the sequences by cloneFraction. 
        Saves the top self.seqs_to_save amount of sequences. 
        Reference:
        Beshnova D, et al. De novo prediction of cancer-associated T cell receptors for noninvasive cancer detection. 
        They implemented 40,000 seqs per class, I am using 50,000 as my default value. 
        '''

        self.logger.info(f"Saving the top {self.seqs_to_save} sequences per repertoire.")

        for cancer_type, cancer_df in data.items():
            # Sort by cloneFraction, take of top self.seqs_to_save sequences. 
            data[cancer_type] = cancer_df.sort_values(by='cloneFraction', ascending=False).iloc[:self.seqs_to_save]
        return data

    # Abstract method.
    def sub_sample_repertoire(self, data:dict) -> dict:
        '''
        Take each repertoire and subsample it self.subsample_size times.
        Returns a dictionary of lists where each key is the disease name and each value is a list of lists.
        Each sub-list contains the sequences sampled from the repertoire.
        '''

        self.logger.info(f"Sub-sampling each repertoire into {self.num_of_samples} samples of size {self.subsample_size} seqs.")
        
        # Collect the top self.seqs_to_save sequences per repertoire.
        top_n_seq_repertoires = self.save_top_n_seqs(data)

        # Iterate through the repertoires and subsample each one.
        repertoire_subsamples = {}
        for repertoire in top_n_seq_repertoires.keys():
            # Dictionary structure: {repertoire: [[subsample1], [subsample2], ...]}
            repertoire_subsamples[repertoire] = [top_n_seq_repertoires[repertoire]['AASeq'].sample(self.subsample_size, replace=True).to_list()
                                                 for i in range(self.num_of_samples)]
            self.logger.debug(f"The number of subsamples per reperotire = {len(repertoire_subsamples[repertoire])}")
            self.logger.debug(f"The number of sequences in each subsample = {len(repertoire_subsamples[repertoire][0])}")

        return repertoire_subsamples
    
    def encode_sequence(self):
        pass
    

    # Abstract method.
    def create_targets(self):
        pass

    # Abstract method.
    def save_data(self):
        self.handler.save_as_csv()

    # Abstract method.
    def complete_preprocess(self):
        repertoires = self.collect_repertoires()
        top_n_seq_repertoires = self.save_top_n_seqs(repertoires)
        self.handler.save_multi_as_csv(top_n_seq_repertoires)