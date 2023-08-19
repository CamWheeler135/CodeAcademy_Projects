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
    """
    KmerPreprocessor performs the pre-processing for the UnsupervisedVisualizer class. 
    
    Preprocessing:
        Collect the data.
        Subsample the repertoire.
        Count the number of k-mers in the subsample.
        Construct a k-mer DataFrame.
        Save the DataFrames. 

    Attributes:
        logger: Logger module for the instance.
        subsample_size: The number of sequences taken from a repertoire in order to conduct k-mer analysis (default value 100 sequences).
        num_of_samples: The number of times 'subsample_size' is taken from a repertoire (default value 1000).
        kmer_size: The length of the k-mer used when conducting k-mer analysis (default value 2).
        base_kmer_counter: Dictionary of the different permutations of k-mer according to k-mer size, this is used to create fresh copies 
                           of the k-mer counter when conducting k-mer analysis.
        handler: DataHandler instance that loads and saves the data into the instantiated paths. 
    """

    def __init__(self, handler: DataHandler, subsample_size:int=100, 
                 num_of_samples:int=1000, kmer_size:int=2):
        
        self.logger = logging.getLogger('Unsupervised Preprocessor')
        self.subsample_size = subsample_size
        self.num_of_samples = num_of_samples
        self.kmer_size = kmer_size
        self.base_kmer_counter = {''.join(kmer): 0 for kmer in product('ARNDCEQGHILKMFPSTWYV', repeat=self.kmer_size)}
        self.handler = handler

    # Abstract method.
    def collect_repertoires(self) -> dict:
        """
        Utilizes the DataHandler that is passed into the preprocessor upon instantiation.
        """
        return self.handler.collect_csv_files()

    # Abstract method.
    def sub_sample_repertoire(self, data:dict[str, pd.DataFrame]) -> dict[str, list[list]]:
        """
        Iterates through the data subsampling 'subsample_size' sequences from the repertoire 'num_of_samples' times.

        Args:
            data: A dictionary of repertoire names (str) and data-frame of sequence information (pd.DataFrame).

        Returns:
            repertoire_subsamples: A dictionary of repertoire names (str) keys with a list value, this list contains sub-lists of sequences sampled 
            from the repertoire with the following structure.

            {'repertoire name 1': [[sequences in here], [sequences in here], [sequences in here]],
            'repertoire name 2': [[sequences in here], [sequences in here], [sequences in here]],
                                                    ...
            'repertoire name x': [[sequences in here], [sequences in here], [sequences in here]]}
        """

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
    def create_targets(self, data:dict[str, pd.DataFrame]) -> dict[str, int]:
        """
        Creates a target value (int) dictionary for the dataset.

        Args: 
            data: A dictionary of repertoire names (str) and data-frame of sequence information (pd.DataFrame).

        Returns: 
            A dictionary where each disease repertoire has an individual numerical target.

            {'repertoire 1 name': 1,
             'repertoire 2 name': 2,
                        ...
             'repertoire x name': x}
        """
        return {disease: value for value, disease in enumerate(sorted(data.keys()))} 
    
    def count_kmers(self, sequences:list[str]) -> dict[str, int]:
        """
        Iterates through each sequence in the sub-sample, counting the number of kmers present.

        Args: 
            sequences: a list of sequence strings [sequence 1, sequence 2, ..., sequence 3]. 

        Returns:
            kmer_vals: A dictionary of k-mers and their counts.

            {'AA': 4, 'AC': 10, ..., 'XX': int}
        """
        # Construct the kmer dictionary.
        kmer_vals = self.base_kmer_counter.copy()
        # Iterate through the sequences in the subsample and increment the dictionary values.
        for sequence in sequences:
            kmer_iterator = nltk.ngrams(sequence, self.kmer_size)
            for kmer in kmer_iterator:
                kmer_vals[''.join(kmer)] += 1 
        return kmer_vals

    def collect_subsample_kmers(self, repertoire_subsamples:dict[str, list[list]]) -> list[dict]:
        """
        Iterate through repertoire subsamples and count the occurrences of k-mers in each subsample.
        For each repertoire, the function will also add the target value (int) at the end of dictionary. 

        For each subsample:
            Count the number of kmers present.
            Add the target value to each repertoire.
            Append the results to a list.
            Return the list.
        
        Args:
            repertoire_subsamples: A dictionary of repertoire name (str) keys with a list of lists where each
            sublist contains sampled sequences. 

        Returns:
            kmer_list: A list of dictionaries of k-mer counts and target values with the following structure:

            [{'AA': 4, 'AC': 10, ..., 'XX': int, 'target': int},
            {'AA': 20, 'AC': 1, ..., 'XX': int, 'target': int},
                                    ...
            {'AA': 4, 'AC': 10, ..., 'XX': int, 'target': int}]

        """

        # Create the targets.
        target_vals = self.create_targets(repertoire_subsamples)
        self.logger.debug(f"Target values: {target_vals}")

        # Count the kmers and add the target to each dict.
        kmer_list = []
        for repertoire in repertoire_subsamples.keys():
            self.logger.info(f"Counting the kmers of size {self.kmer_size} present in each {repertoire} sub-sample.'")
            for subsample in repertoire_subsamples[repertoire]:
                subsample_kmers = self.count_kmers(subsample) # Count the occurrences of k-mers in the subsample. 
                # Adds the target value to the dictionary of kmers.
                subsample_kmers['target'] = target_vals[repertoire]
                kmer_list.append(subsample_kmers)

        return kmer_list

    # Abstract method.
    def construct_kmer_df(self, kmer_list:list[dict]) -> pd.DataFrame:
        """
        Forms a pd.DataFrame out of all the repertoires k-mer dictionaries.
        
        Args:
            kmer_list: A list of dictionaries where each dictionary is a count of the k-mers present in
                       a subsample and its target integer value.

        Returns:
            kmer_df: A pandas DataFrame, rows correspond to a sample, columns correspond to a k-mer.
        """
        kmer_df = pd.DataFrame(kmer_list)
        self.logger.debug(f"The shape of kmer data frame is {kmer_df.shape}")
        return kmer_df
    
    def save_data(self, kmer_df: pd.DataFrame):
        """ 
        Save the data frame as a CSV file in the handlers 'save_path' directory.
        """
        self.handler.save_as_csv(kmer_df)

    # Abstract method.
    def complete_preprocess(self):
        """
        Preprocess the data into a format that can be used by the UnsupervisedVisualizer class.
        1. Load in the repertoires.
        2. Subsample each repertoire.
        3. Create the desired k-mers.
        4. Take the sub-samples and count the kmers present in each subsample.
        5. Construct a pd.DataFrame of the k-mer counts in each subsample.
        6. Save the data as a CSV.
        """
        repertoires = self.collect_repertoires()
        subsamples = self.sub_sample_repertoire(repertoires)
        kmers = self.collect_subsample_kmers(subsamples)
        kmer_df = self.construct_kmer_df(kmers)
        self.save_data(kmer_df)


class CNNPreprocessor(BasePreprocessor):
    """
    CNNPreprocessor class performs all of the preprocessing needed for training the CNN model.
    Note: This class is NOT responsible for loading the data for the model, only to take the raw sequence data
          and preprocess it into biophysiochemical features. The DataHandler that is passed into the CNNPreprocessor
          upon instantiation is then responsible for saving the data into the appropriate directory structure.

    Preprocessing:
    
    Attributes:
        logger: Logger for the instance.
        handler: DataHandler instance that loads and saves the data into the instantiated paths. 
        seqs_to_save: The value controls the top number of sequences that are saved from each repertoire (based on CloneFraction).
        subsample_size: The number of sequences that are included in a subsample.
        num_of_samples: The number of subsamples taken from each repertoire.
        encoding_matrix: The matrix that contains the biophysiochemcial features for each amino acid, this is used to
                         encode the sequence data for the CNN.
    """

    def __init__(self, handler: DataHandler, seqs_to_save:int=50000, 
                 subsample_size:int=100, num_of_samples:int=40000):
        
        self.logger = logging.getLogger('CNN Preprocessor')
        self.handler = handler
        self.seqs_to_save = seqs_to_save
        self.num_of_samples = num_of_samples
        self.subsample_size = subsample_size
        self.encoding_matrix = None # Add the Beshova matrix here. 

    def collect_repertoires(self):
        """
        Utilizes the DataHandlers method to load the data from csv to DataFrames.
        """
        return self.handler.collect_csv_files()

    def save_top_n_seqs(self, data:dict[str, pd.DataFrame]) -> dict:
        """
        Takes each repertoire and sorts the sequences by cloneFraction before
        saving the top 'seqs_to_save' amount of sequences. 
        Reference:
        Beshnova D, et al. De novo prediction of cancer-associated T cell receptors for noninvasive cancer detection. 
        They implemented 40,000 seqs per class, I am using 50,000 as my default value. 

        Args:
            data: A dictionary of repertoire names (str) and data-frame of sequence information (pd.DataFrame).
        
        Returns:
            data: The exact same structure as the argument, but contains only the top n sequences.

        """

        self.logger.info(f"Saving the top {self.seqs_to_save} sequences per repertoire.")

        for cancer_type, cancer_df in data.items():
            # Sort by cloneFraction, take off top self.seqs_to_save sequences. 
            data[cancer_type] = cancer_df.sort_values(by='cloneFraction', ascending=False).iloc[:self.seqs_to_save]
        return data

    def sub_sample_repertoire(self, data:dict[str, pd.DataFrame]) -> dict[str, list[list]]:
        """
        Iterates through the data subsampling 'subsample_size' sequences from the repertoire 'num_of_samples' times.

        Args:
            data: A dictionary of repertoire names (str) and data-frame of sequence information (pd.DataFrame).

        Returns:
            repertoire_subsamples: A dictionary of repertoire names (str) keys with a list value, this list contains sub-lists of sequences sampled 
            from the repertoire with the following structure.

            {'repertoire name 1': [[sequences in here], [sequences in here], [sequences in here]],
            'repertoire name 2': [[sequences in here], [sequences in here], [sequences in here]],
                                                    ...
            'repertoire name x': [[sequences in here], [sequences in here], [sequences in here]]}
        """

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