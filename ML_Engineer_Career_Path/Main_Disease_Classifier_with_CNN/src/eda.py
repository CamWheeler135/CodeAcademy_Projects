''' File that contains the code for exploratory data analysis. '''

# Module Imports
from pathlib import Path
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

# Local Imports
from src.base_logger import logging
from src.data_cleaning import DataHandler

class ExploratoryAnalyzer():
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
    
    def __init__(self, handler:DataHandler):
        self.logger = logging.getLogger("EDA Analyzer")
        self.handler = handler
        self.sequence_counts = {}
        self.seq_len_metrics = {}
        self.gene_usage = {}
        self.jaccard_index = {}
        self.aa_counts = {}

    def __str__(self):
        return "Exploratory Data Analysis Object."
    
    def collect_repertoires(self) -> dict:
         return self.handler.collect_csv_files()
    
    # Collect the number of sequences we have for each disease.
    def count_sequences_per_disease(self, data:dict) -> dict:
        ''' Counts the number of sequences present in each repertoire. '''

        self.logger.info("Counting the number of sequences in each repertoire.")
        for repertoire in data:
            self.sequence_counts[repertoire] = len(data[repertoire].index)
        
    def count_sequence_length_metrics(self, data:dict) -> dict:
        ''' 
        Computes metrics regarding the length of sequences in the repertoires.
        - Seq len distribution (repertoire).
        - Mean seq len (dataset).
        - Std seq len (dataset).
        '''

        self.logger.info("Computing the sequence length metrics.")

        # Variables for computing database mean length and standard deviation.
        total_seq_lens = 0
        total_num_of_seqs = 0
        repertoire_seq_lens = []

        # Computing the sequence lengths for each repertoire.
        for repertoire in data:
            lengths = data[repertoire]['AASeq'].str.len()
            self.seq_len_metrics[repertoire] = {'len_counts': lengths.value_counts(normalize=True)}
            # Manipulate variables for computing database mean length and standard deviation.
            repertoire_seq_lens.append(list(lengths)) 
            total_seq_lens += np.sum(lengths)
            total_num_of_seqs += self.sequence_counts[repertoire]
            self.logger.debug(self.seq_len_metrics[repertoire])

        # Computing the mean.
        self.logger.debug(f"The total length of sequences in the database are {total_seq_lens}")
        self.logger.debug(f"The overall number of sequences in the database are {total_num_of_seqs}")
        mean = total_seq_lens / total_num_of_seqs
        self.seq_len_metrics['mean_len'] = mean

        # Computing the standard deviation. 
        len_array = np.hstack(repertoire_seq_lens)
        standard_dev = np.sqrt((np.sum((len_array - mean) ** 2)) / total_num_of_seqs)
        self.seq_len_metrics['standard_dev'] = standard_dev
        self.logger.debug(f"The standard deviation is {standard_dev}")

    def count_gene_usage(self, data:dict):
        ''' Counts each permutation of gene present in V, D and J region. '''

        gene_regions = ['Vregion', 'Jregion' ,'Dregion']

        # For each repertoire, count the number of permutations in each gene region. 
        for gene in gene_regions:
            self.gene_usage[gene] = {repertoire: data[repertoire][gene].value_counts(normalize=True) for repertoire in data}
            self.logger.debug(f" The usage counts for {gene} are {self.gene_usage[gene]}")
    
    def compute_jaccard_index(self, data:dict):
        ''' Computes the Jaccard similarity index between each of the diseases. '''

        self.logger.info("Computing the Jaccard Index.")
        for repertoire_A in data:    

            # Holds the jaccard score for repertoire A vs repertoire B
            jaccard_scores = {}

            for repertoire_B in data:

                # No point in comparing it to itself. 
                if repertoire_A == repertoire_B:
                    continue
                # Compute the Jaccard index between repertoire A and all other repertoires in the dataset. 
                intersection = len(set(data[repertoire_A]['AASeq']).intersection(data[repertoire_B]['AASeq']))
                self.logger.debug(f"The intersection between {repertoire_A} and {repertoire_B} is {intersection}.")
                union = len(set(data[repertoire_A]['AASeq']).union(data[repertoire_B]['AASeq']))
                self.logger.debug(f"The union between {repertoire_A} and {repertoire_B} is {union}.")
                jaccard_scores[repertoire_B] = intersection/union
            
            # Dictionary hold the jaccard index between repertoire A and all other repertoires. 
            self.jaccard_index[repertoire_A] = jaccard_scores
            self.logger.debug(f"Jaccard index for {repertoire_A} = {jaccard_scores}")

    def count_aa_occurrences(self, data:dict):
        ''' Counts the proportion of occurrences each amino acid has in a repertoire. '''
        self.logger.info("Counting amino acids present in the repertoire.")
        for repertoire in data:
            # Splits the sequence into characters making a new dataframe.
            # First column and '' at end of seq needs to be removed due to quirk in str.split method in pandas. 
            aa_breakdown = (data[repertoire]['AASeq'].str.split('', expand=True).iloc[:, 1:-1].replace({'': None}))
            self.aa_counts[repertoire] = {position : aa_breakdown[position].value_counts(normalize=True) for position in aa_breakdown.columns}
            self.logger.debug(self.aa_counts[repertoire])
    
    def complete_eda(self):
        ''' Performs the complete EDA. '''
        repertoires = self.collect_repertoires()
        self.count_sequences_per_disease(repertoires)
        self.count_sequence_length_metrics(repertoires)
        self.count_gene_usage(repertoires)
        self.compute_jaccard_index(repertoires)
        self.count_aa_occurrences(repertoires)

class UnsupervisedVisualizer():
    '''
    Class will perform unsupervised learning for visualization.
    Data preprocessing is performed by the Unsupervised Preprocessor class.
    Utilizes Sklearn's pipeline class to perform analysis.
    Performs PCA analysis projecting the data into 2D.
    Plots are created using the plotter object. 
    '''

    def __init__(self, handler):
        self.logger = logging.getLogger("Unsupervised Visualizer.")
        self.handler = handler
        self.pca_data = None
        self.pca_targets = None

    def collect_kmer_dataset(self):
        return self.handler.collect_csv_files()

    def perform_pca(self, data:dict):
        ''' Performs PCA analysis on the data. '''

        self.logger.info("Performing PCA Analysis.")

        pca_pipeline = Pipeline([("scaler", StandardScaler()),
                                ("pca", PCA(n_components=2))])
        
        assert len(data.keys()) == 1 # Only one dataset should be present.
        
        for kmer_dataset in data.values():
            X_train = kmer_dataset.drop(columns='target')
            Y_train = kmer_dataset['target']
            self.pca_data = pd.DataFrame(pca_pipeline.fit_transform(X_train))
            self.pca_targets = Y_train
        
    def save_pca_df(self):
        self.handler.save_as_csv(self.pca_data)

    def complete_pca_analysis(self):
        repertoires = self.collect_kmer_dataset()
        self.perform_pca(repertoires)
        self.save_pca_df()


class Plotter():
    '''
    Class interacts with ExploratoryAnalyzer and UnsupervisedVisualizer objects to plot visualizations. 
    
    Plots:
        - Bar plot of sequence counts per cancer type.
        - Pie chart of sequences per cancer type in the whole dataset.
        - Lineplot of sequence length distribution per cancer type.
        - Heatmap of gene usage per cancer type.
        - Heatmap of Jaccard index between cancer types.
        - Heatmap of proportional amino acid presence per cancer type.
    '''

    def __init__(self,eda: ExploratoryAnalyzer, unsupervised_vis: UnsupervisedVisualizer):
        self.logger = logging.getLogger("Plotter")
        self.eda = eda
        self.unsupervised_vis = unsupervised_vis
        self.eda_save_path = 'Figures/EDA/'
        self.unsupervised_save_path = 'Figures/Unsupervised/'
    
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
        plt.savefig(Path(self.eda_save_path + 'Sequence_Counts'))
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
        plt.savefig(Path(self.eda_save_path + 'Cancer_Type_Proportions'))
        self.logger.info("Database proportion pie chart saved.")

    def plot_sequence_len_distribution(self):
        '''Plots the sequence length distribution of each repertoire.'''

        # Data.
        mean = self.eda.seq_len_metrics['mean_len']
        standard_dev = self.eda.seq_len_metrics['standard_dev']
        del self.eda.seq_len_metrics['mean_len']
        del self.eda.seq_len_metrics['standard_dev']

        # Plot.
        fig = plt.figure(figsize=(10, 10))
        for repertoire in self.eda.seq_len_metrics.keys():
            sns.lineplot(x=self.eda.seq_len_metrics[repertoire]['len_counts'].index, 
                         y=self.eda.seq_len_metrics[repertoire]['len_counts'], label=str(repertoire))
        plt.axvline(mean , color='k', ls='--', lw=2)
        plt.axvline(mean - standard_dev, color='k', ls='--', lw=1)
        plt.axvline(mean + standard_dev, color='k', ls='--', lw=1)
        plt.title("Proportional Sequence Length")
        plt.ylabel("Proportion of Repertoire")
        plt.xlabel("Length of Sequence")
        plt.legend()
        plt.savefig(Path(self.eda_save_path + 'Amino_Acid_Distributions'))
        self.logger.info("Sequence length distribution chart saved.")

    def plot_gene_usage(self):
        ''' Plots a heatmap of the gene permutations in each disease for V, D and J region '''

        # There are only 3 regions in analysis so hard coding the number of graphs should be okay. 
        # Collect the data for each of the heatmaps.
        vgene_df = pd.DataFrame(self.eda.gene_usage['Vregion']).fillna(0).transpose()
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
        plt.savefig(Path(self.eda_save_path + 'Gene Usage Heatmaps'))

    def plot_jaccard_index(self):

        # Data.
        repertoires = self.eda.jaccard_index.keys()
        jaccard_value_list = []
        for repertoire in repertoires:
            jaccard_value_list.append(self.eda.jaccard_index[repertoire])
        jaccard_df = pd.DataFrame(jaccard_value_list, index=repertoires, columns=repertoires)
        self.logger.debug(jaccard_df)

        # Plot.
        fig = plt.figure(figsize=(10, 10))

        mask = np.zeros(jaccard_df.shape, dtype=bool) 
        mask[np.triu_indices(len(mask))] = True 
        np.fill_diagonal(jaccard_df.values, 0)
        self.logger.debug(jaccard_df)

        sns.heatmap(jaccard_df, vmin=0.001, vmax=0.2, mask=mask)
        plt.tight_layout()
        plt.savefig(Path(self.eda_save_path + 'Jaccard_Heatmap'))


    # Plot the amino acid distribution for each repertoire (heatmap).
    def plot_aa_distribtuion(self):
        ''' Plots the amino acid distribution for each repertoire in a heatmap. '''

        # Plot.
        keys = list(self.eda.aa_counts.keys())
        fig, ax = plt.subplots(3, 3, figsize=(30, 20), sharey='row')
        cbar_ax = fig.add_axes([.91, .3, .03, .4])

        for i in range(len(keys)):
            data = pd.DataFrame(self.eda.aa_counts[keys[i]]).fillna(0)
            sns.heatmap(data, vmax=0.4, cbar=i == 0, ax=ax.flat[i], cbar_ax=cbar_ax)
            ax.flat[i].set_title(f"{keys[i]} Amino Acid Distribution")
        fig.tight_layout(rect=[0, 0, .9, 1])
        plt.savefig(Path(self.eda_save_path + 'Amino_Acid_Distributions'))
    
    def plot_pca(self):
        ''' Plots a scatter graph of the pca data. '''

        # Data.
        pca_x= self.unsupervised_vis.pca_data.iloc[:, 0]
        pca_y = self.unsupervised_vis.pca_data.iloc[:, 1]
        pca_targets = self.unsupervised_vis.pca_targets

        # Plot.
        fig = plt.figure(figsize=(10, 10))
        plt.scatter(pca_x, pca_y, c=pca_targets, cmap='tab10')
        plt.title("PCA Scatter Plot")
        plt.xlabel("PCA 1")
        plt.ylabel("PCA 2")
        plt.savefig(Path(self.unsupervised_save_path + 'PCA_Scatter_Plot'))

    def complete_plots(self):
        ''' Runs all EDA and Unsupervised plots. '''

        self.plot_sequence_count_bar()
        self.plot_sequence_count_pie()
        self.plot_sequence_len_distribution()
        self.plot_gene_usage()
        self.plot_jaccard_index()
        self.plot_aa_distribtuion()
        self.plot_pca()
        self.logger.info("All plots saved.")