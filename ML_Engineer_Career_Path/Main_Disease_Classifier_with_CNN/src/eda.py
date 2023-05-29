''' File that contains the code for exploratory data analysis. '''

# Module Imports
from pathlib import Path
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Local Imports
from src.base_logger import logging

class DataLoader():
    '''
    Class will load the data into a dictionary for use in analysis.
    '''

    # Function that loads each of the CSV files for analysis.
    def __init__(self, collection_path):
        self.collection_path = collection_path

    def collect_files(self):
        # Returns the file name (stripping the csv tag) with a data frame of the sequences.
        return {os.path.splitext(file.name)[0]: pd.read_csv(file) for file in self.collection_path.iterdir()}

class ExploratoryAnalyzer(DataLoader):
    '''
    Class will extract and store information for exploratory data analysis. 
    All plots are created using the 'Plotter" object. 

    Data Analysis:
        - Number of sequences per disease.
        - Sequence length. 
        - Sequence sharing (Jacquard Index).
        - Gene usage. 
        - Repertoire amino acid distributions.
    '''
    
    def __init__(self, collection_path=Path('data/cleaned_files')):
        super().__init__(collection_path)
        self.logger = logging.getLogger("EDA Analyzer")
        self.sequence_counts = {}
        self.seq_len_info = {}
        self.gene_usage = {}

    def __str__(self):
        return "Exploratory Data Analysis Object."
    
    # Collect the number of sequences we have for each disease.
    def count_sequences_per_disease(self, data:dict) -> dict:
        ''' Counts the number of sequences present in each repertoire. '''

        self.logger.info("Counting the number of sequences in each repertoire.")
        for repertoire in data:
            self.sequence_counts[repertoire] = len(data[repertoire].index)
        
    def count_sequence_length_mertics(self, data:dict) -> dict:
        ''' Computes metrics regarding the length of sequences in the repertoires.'''

        self.logger.info("Computing the sequence length metrics.")

        # Variables for computing database mean length and standard deviation.
        total_seq_lens = 0
        total_num_of_seqs = 0
        repertoire_seq_lens = []

        # Computing the sequence lengths for each repertoire.
        for repertoire in data:
            lengths = data[repertoire]['AASeq'].str.len()
            repertoire_seq_lens.append(list(lengths)) 
            self.seq_len_info[repertoire] = {'len_counts': lengths.value_counts() / len(data[repertoire].index) * 100}
            total_seq_lens += np.sum(lengths)
            total_num_of_seqs += self.sequence_counts[repertoire]
            self.logger.debug(self.seq_len_info[repertoire])

        # Computing the mean.
        self.logger.debug(f"The total length of sequences in the database are {total_seq_lens}")
        self.logger.debug(f"The overall number of sequences in the database are {total_num_of_seqs}")
        mean = total_seq_lens / total_num_of_seqs
        self.seq_len_info['mean_len'] = mean

        # Computing the standard deviation. 
        len_array = np.hstack(repertoire_seq_lens)
        standard_dev = np.sqrt((np.sum((len_array - mean) ** 2)) / total_num_of_seqs)
        self.seq_len_info['standard_dev'] = standard_dev
        self.logger.debug(f"The standard deviation is {standard_dev}")

    def count_gene_usage(self, data:dict):
        ''' Counts each permutation of gene present in V, D and J region. '''

        gene_regions = ['Vregion', 'Jregion' ,'Dregion']

        # For each repertoire, count the number of permutations in each gene region. 
        for gene in gene_regions:
            self.gene_usage[gene] = {repertoire: data[repertoire][gene].value_counts(normalize=True) for repertoire in data}
            self.logger.debug(f" The usage counts for {gene} are {self.gene_usage[gene]}")

        '''TODO'''
        # Comparing the number of unknowns in the repertoire. 

    
    '''TODO'''
    # Collect Jaccard Index of each disease.

    # Collect the gene usage of each disease. 
        # Again we want a proportion here, what proportion sequences in a repertoire use this gene. 

    # Collect the Amino acid distribution for each disease. 
        # We want a proportion here, otherwise we will just be getting high scores for large repertoires. 



class UnsupervisedVisualizer(DataLoader):
    '''
    Class will perform all of the preprocessing required to perform unsupervised learning.
    Utilizes a 'Bag of Words' approach, counting the occurrences of k-mer in each subsample.
    Creates a data frame of the counts per sample then performs PCA analysis projecting the data into 2D.
    Plots are created using the plotter object. 
    '''

    def __init__(self):
        self.logger = logging.getLogger("Unsupervised Visualizer.")

    # Perform Unsupervised learning.
        # Subsample each file. 
        # Count the Kmers in each subsample
        # PCA into 2 dimensions, be sure to include colors so we can visually see what is happening. 


class Plotter():
    '''
    Class interacts with ExploratoryAnalyzer object to plot visualizations. 
    
    Plots:
        - 
    '''

    def __init__(self, eda, save_path='Figures/EDA/'): # unsupervised_vis
        self.logger = logging.getLogger("Plotter")
        self.eda = eda
        # self.unsupervised_vis = unsupervised_vis
        self.save_path = save_path
    
    def plot_sequence_count_bar(self):
        ''' Plots the number of sequences each cancer type contains. '''

        fig = plt.figure(figsize=(12,15))
        sns.barplot(x=list(self.eda.sequence_counts.keys()), y=list(self.eda.sequence_counts.values()))
        plt.xlabel("Cancer Types")
        plt.ylabel("Number of Sequences")
        plt.xticks(fontsize=10, rotation=45)
        plt.yticks(fontsize=10)
        plt.ticklabel_format(style='plain', axis='y')
        plt.title("Comparison of Total Sequence Count per Cancer Type")
        plt.savefig(Path(self.save_path + 'Sequence_Counts'))
        self.logger.info("Sequence bar chart saved.")

    def plot_sequence_count_pie(self):
        '''Plots the proportion of each cancer type in the database'''

        # Data manipulation. 
        repertoire_seqs = list(self.eda.sequence_counts.values())
        total_seqs = np.sum(repertoire_seqs)
        self.logger.debug(f"Total number of sequences in the database is {total_seqs}")
        self.logger.debug(f"The values of each repertoire are {repertoire_seqs}")
        percent_of_database = [(value / total_seqs) * 100 for value in repertoire_seqs]
        self.logger.debug(f"The percentages of each cancer type are {percent_of_database}")

        # Plotting code
        fig = plt.figure(figsize=(10,10))
        plt.pie(percent_of_database, labels=list(self.eda.sequence_counts.keys()), autopct='%1.1f%%')
        plt.title("Percentage of Cancer Type in Total Database.")
        plt.savefig(Path(self.save_path + 'Cancer_Type_Proportions'))
        self.logger.info("Database proportion pie chart saved.")

    def plot_sequence_len_distribution(self):
        '''Plots the sequence length distribution of each repertoire.'''

        # Data.
        mean = self.eda.seq_len_info['mean_len']
        standard_dev = self.eda.seq_len_info['standard_dev']
        del self.eda.seq_len_info['mean_len']
        del self.eda.seq_len_info['standard_dev']

        # Plot.
        fig = plt.figure(figsize=(10, 10))
        for repertoire in self.eda.seq_len_info.keys():
            sns.lineplot(x=self.eda.seq_len_info[repertoire]['len_counts'].index, y=self.eda.seq_len_info[repertoire]['len_counts'], label=str(repertoire))
        plt.axvline(mean , color='k', ls='--', lw=2)
        plt.axvline(mean - standard_dev, color='k', ls='--', lw=1)
        plt.axvline(mean + standard_dev, color='k', ls='--', lw=1)
        plt.title("Proportional Sequence Length")
        plt.ylabel("Proportion of Repertoire")
        plt.xlabel("Length of Sequence")
        plt.legend()
        plt.savefig(Path(self.save_path + 'Amino_Acid_Distributions'))
        self.logger.info("Sequence length distribution chart saved.")


    def plot_gene_usage(self):
        ''' Plots a heatmap of the gene permutations in each disease for V, D and J region '''

        # There are only 3 regions in analysis so hard coding the number of graphs should be okay. 
        # Collect the data for each of the heatmaps.
        vgene_df = pd.DataFrame(self.eda.gene_usage['Vregion']).dropna(axis=0, how='all').fillna(0).transpose()
        dgene_df = pd.DataFrame(self.eda.gene_usage['Dregion']).fillna(0).transpose()
        jgene_df = pd.DataFrame(self.eda.gene_usage['Jregion']).fillna(0).transpose()
        self.logger.debug(vgene_df)
        self.logger.debug(dgene_df)
        self.logger.debug(jgene_df)
        
        # Plots
        fig, ax = plt.subplots(3, 1, figsize=(30, 25), constrained_layout=True)
        ax[0] = sns.heatmap(vgene_df, ax=ax[0])
        ax[0].set_title("Percentage V Gene Usage Between Diseases.")
        ax[1] = sns.heatmap(dgene_df, ax=ax[1])
        ax[1].set_title("Percentage D Gene Usage Between Diseases.")
        ax[2] = sns.heatmap(jgene_df, ax=ax[2])
        ax[2].set_title("Percentage J Gene Usage Between Diseases.")
        plt.xticks(fontsize=8, rotation=90)
        plt.yticks(fontsize=15)
        plt.savefig(Path(self.save_path + 'Gene Usage Heatmaps'))



    # Plot the Jacquard Index (Heatmap).

    # Plot the amino acid distribution for each repertoire (Sequence logo graph).

    # Plot the unsupervised learning charts, this should allow us to see if subsample's from the different cancer types cluster together. 
        # Scatter graph, remember to color them according to disease so we can see what is going where. 

