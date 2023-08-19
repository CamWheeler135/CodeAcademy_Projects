''' File that contains the code for exploratory data analysis. '''

# Module Imports
from pathlib import Path
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

# Local Imports
from src.base_logger import logging
from src.data_cleaning import DataHandler

class ExploratoryAnalyzer():
    """
    The ExploratoryAnalyzer will perform analysis on CDR3β dataset and store the extracted information as attributes.
    All analysis plots are created using a 'Plotter" object that accesses these attributes. 

    Analysis:
        Number of sequences present in each repertoire.
        Sequence length metrics of each repertoire.
        V, D and J region gene usage of each repertoire.
        Jaccard index comparing sequence similarity between repertoires.
        Amino acid positioning in each repertoire.

    Attributes:
        logger: Logger instance for the class.
        handler: DataHandler instance that handles collecting data from given directory.
        sequence_counts: A dictionary that contains the number of sequences found in each repertoire.
        seq_len_metrics: A dictionary that contains the number of times a sequence of length 'x' appears in the repertoire,
                         At the end of the dictionary there are dataset metrics -> mean sequnce length and 
                         sequence length standard distribution.
        gene_usage: A dictionary that contains the number of times a certain gene type is found in the V region, D region 
                    and J region of the CDR3 sequence.
        aa_counts: A dictionary that contains the proportion of times an amino acid appears in at each position 1 -> n of a sequence
                   for each repertoire.
    """
    
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
        """ 
        Collects repertoire CSV files using the DataHandler that is passed into the EDA
        object as an argument.
        """
        return self.handler.collect_csv_files()
    
    def count_sequences_per_disease(self, data:dict[str, pd.DataFrame]) -> None:
        """ 
        Counts the number of sequences present in each repertoire. 
        
            Args:
            data: A dictionary of repertoire names (str) and data-frame of sequence information (pd.DataFrame).
            
            Returns:
                None

            Sets the sequence_count attribute of the instance to a dictionary of repertoire names (str) keys 
            and number of sequences in the repertoire (int) values. 
                {'repertoire 1 name': 10,
                 'repertoire 2 name': 44,
                            ...
                 'repertoire x name: 112000}
        """

        self.logger.info("Counting the number of sequences in each repertoire.")
        for repertoire in data:
            seq_count = len(data[repertoire].index)
            self.sequence_counts[repertoire] = seq_count
            self.logger.debug(f"Number of sequences in repertoire {repertoire} ----> {seq_count}")
        
    def count_sequence_length_metrics(self, data:dict[str, pd.DataFrame]) -> None:
        """
        Computes metrics regarding the length of sequences in the repertoires.
        1. Counts the occurrences of each sequence length.
        2. Dataset mean sequence length.
        3. Dataset sequence length standard deviation.
        Sets the value of the instances seq_len_metrics attribute.

        Args:
            data: A dictionary of repertoire names (str) and data-frame of sequence information (pd.DataFrame).

        Returns:
            None

        Sets the seq_len_metrics attribute of the instance to a dictionary of repertoire sequence length counts, 
        alongside the dataset mean and standard deviation with the following structure:

            {'repertoire 1 name': {'len_counts': },
             'repertoire 2 name': {'len_counts: },
                            ...
             'repertoire x name': {'len_counts: }.
             'mean_len': 14.222,
             'standard_dev': 1.777}
        """

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
        self.logger.debug(f"The total length of sequences in the database ----> {total_seq_lens}")
        self.logger.debug(f"The overall number of sequences in the database ----> {total_num_of_seqs}")
        mean = total_seq_lens / total_num_of_seqs
        self.seq_len_metrics['mean_len'] = mean

        # Computing the standard deviation. 
        len_array = np.hstack(repertoire_seq_lens)
        standard_dev = np.sqrt((np.sum((len_array - mean) ** 2)) / total_num_of_seqs)
        self.seq_len_metrics['standard_dev'] = standard_dev
        self.logger.debug(f"The datasets standard deviation ----> {standard_dev}")

    def count_gene_usage(self, data:dict[str, pd.DataFrame]) -> None:
        """ 
        Counts the proportion of each permutation of gene present in V, D and J region of each sequence in each repertoire. 
        Sets the value of the instances gene_usage attribute.

        Args:
            data: A dictionary of repertoire names (str) and data-frame of sequence information (pd.DataFrame).

        Returns: 
            None

        Sets the gene_usage attribute of the instance to a dictionary of region names (str) keys with a dictionary as the value 
        containing repertoire name (str) keys and gene usage proportion counts (pd.Series) values with the following structure:

        {'Vregion': {'repertoire 1 name': pd.Series,
                      'repertoire 2 name': pd.Series,
                                ...
                      'repertoire x name': pd.Series},

        'Dregion': {'repertoire 1 name': pd.Series,
                     'repertoire 2 name': pd.Series,
                                ...
                     'repertoire x name': pd.Series},

        'Jregion': {'repertoire 1 name': pd.Series,
                     'repertoire 2 name': pd.Series,
                                ...
                     'repertoire x name': pd.Series}}
        """

        gene_regions = ['Vregion', 'Jregion' ,'Dregion']

        # For each repertoire, count the number of permutations in each gene region. 
        for gene in gene_regions:
            self.gene_usage[gene] = {repertoire: data[repertoire][gene].value_counts(normalize=True) for repertoire in data}
            self.logger.debug(f" The usage counts for {gene} are {self.gene_usage[gene]}")
    
    def compute_jaccard_index(self, data:dict[str, pd.DataFrame]) -> None:
        """ 
        Computes the Jaccard similarity index between each of the disease repertoires in the dataset. 
        Sets the value of the instances jaccard_index attribute.

        Args:
            data: A dictionary of repertoire names (str) and data-frame of sequence information (pd.DataFrame).
        
        Returns: 
            None

        Sets the 'jaccard_index' attribute of the instance to a dictionary of repertoire A name (str) keys with a dictionary
        as the value, this dictionary contains repertoire B names (str) keys and jaccard score (float) values between 
        repertoire A and B with the following structure:

        {'repertoire 1 name': {'repertoire 2 name': 0.32,
                               'repertoire 3 name': 0.82,
                                        ...
                               'repertoire x name': 0.04},

         'repertoire 2 name': {'repertoire 1 name': 0.43,
                               'repertoire 3 name': 0.34,
                                        ...
                               'repertoire x name': 0.60},

         'repertoire x name': {'repertoire 1 name': 0.50,
                               'repertoire 2 name': 0.30,
                                        ...
                               'repertoire y name': 0.009}}
        """

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
            self.logger.debug(self.jaccard_index)

    def count_aa_occurrences(self, data:dict[str, pd.DataFrame]) -> None:
        """ 
        Counts the proportion of occurrences each amino acid has in a repertoire. 

        Args:
            data: A dictionary of repertoire names (str) and data-frame of sequence information (pd.DataFrame).
        
        Returns:
            None

        Sets the 'aa_counts' attribute of the instance to a dictionary of repertoire name (str) keys with 
        a dictionary as a value containing position in the sequence (int) and the proportion each amino acid occurs 
        in that position (pd.Series) with the following structure:

        {'repertoire 1 name': {position 1: pd.Series,
                               position 2: pd.Series,
                               position 3: pd.Series,
                                        ...
                               position x: pd.Series},

         'repertoire 2 name': {position 1: pd.Series,
                               position 2: pd.Series,
                               position 3: pd.Series,
                                        ...
                               position x: pd.Series},

         'repertoire 3 name': {position 1: pd.Series,
                               position 2: pd.Series,
                               position 3: pd.Series,
                                        ...
                               position x: pd.Series}}
        """

        self.logger.info("Counting amino acids present in the repertoire.")
        for repertoire in data:
            '''TODO''' # This aa_breakdown is pretty much unreadable, need to refactor.
            aa_breakdown = (data[repertoire]['AASeq'].str.split('', expand=True).iloc[:, 1:-1].replace({'': None}))
            self.aa_counts[repertoire] = {position : aa_breakdown[position].value_counts(normalize=True) for position in aa_breakdown.columns}
        self.logger.debug(self.aa_counts)
    
    def complete_eda(self) -> None:
        """ 
        Performs the complete EDA using all of the methods available to the EDA instance.
        1. Collects the data using the DataHandler object passed in at instantiation.
        2. Counts the number of sequences in each repertoire.
        3. Counts the occurances of each sequence length in the repertoire.
        4. Counts the proportion of gene usage in each repertoire.
        5. Computes the Jaccard index between all repertoires.
        6. Counts the occurrences of each amino acid at each position in the sequence for all sequences in the repertoire. 

        Args: 
            None
        
        Returns: 
            None
        """
        repertoires = self.collect_repertoires()
        self.count_sequences_per_disease(repertoires)
        self.count_sequence_length_metrics(repertoires)
        self.count_gene_usage(repertoires)
        self.compute_jaccard_index(repertoires)
        self.count_aa_occurrences(repertoires)


class UnsupervisedVisualizer():
    '''
    UnsupervisedVisualizer will perform Principle Component Analysis unsupervised learning for visualization of repertoire
    amino-acid k-mer data. Utilizes Sklearn's pipeline class to perform analysis, performing PCA analysis projecting 
    the k-mer data into 'num_of_components' dimension (default 2D).

    Data preprocessing MUST be performed by the 'KmerPreprocessor' class prior to running the UnsupervisedVisualizer.
    Plots are created using the 'Plotter' object. 

    Attributes:
        logger: Logger instance for the class.
        handler: DataHandler instance that handles collecting data from given directory.
        num_of_components: The number of components the PCA will reduce too (recommended 2D or 3D for plotting).
        pca_data: Set to none upon instantiation, will be replaced by a pd.DataFrame of PCA data once PCA has been performed.
        pca_targets: The targets associated with the PCA instances. 

    '''

    def __init__(self, handler:DataHandler, num_of_components:int=2):
        self.logger = logging.getLogger("Unsupervised Visualizer.")
        self.handler = handler
        self.num_of_components = num_of_components
        self.pca_data = None
        self.pca_targets = None

    def collect_kmer_dataset(self):
        """
        Collects k-mer CSV files using the DataHandler that is passed into the 
        UnsupervisedVisualizer object as an argument. 
        """

        return self.handler.collect_single_csv_file()

    def perform_pca(self, kmer_dataset:pd.DataFrame):
        """
        Performs PCA analysis on the data.
        
        Args:
            kmer_dataset: A data-frame of k-mer data that is produced prior to analysis by the KmerPreprocessor object.

        Returns:
            None

        Sets the instances pca_data and pca_targets as a pd.DataFrame of reduced_dimension k-mer data 
        and their associated targets respectively.
        """

        self.logger.info("Performing PCA Analysis.")

        pca_pipeline = Pipeline([("scaler", StandardScaler()),
                                ("pca", PCA(n_components=self.num_of_components))])
        
        X_train = kmer_dataset.drop(columns='target')
        Y_train = kmer_dataset['target']
        self.pca_data = pd.DataFrame(pca_pipeline.fit_transform(X_train))
        self.pca_targets = Y_train
        
    def save_pca_df(self):
        """
        Uses the DataHandler object method to save
        the reduced dimensionality k-mer data.
        """
        self.handler.save_as_csv(self.pca_data)

    def complete_pca_analysis(self):
        """
        Performs complete PCA analysis.
        1. Collects the data that has been generated by 'KmerPreprocessor' object.
        2. Performs PCA analysis, producing a pd.DataFrame.
        3. Saves the pd.DataFrame as a CSV file.
        """
        kmer_df = self.collect_kmer_dataset()
        self.perform_pca(kmer_df)
        self.save_pca_df()


class Plotter():
    '''
    Plotter class interacts with ExploratoryAnalyzer and UnsupervisedVisualizer objects to plot visualizations.
    This was designed with singularity in mind, instead of the EDA and UnsupervisedVisualizer objects doing their plots,
    they are instead passed to this class that should extract the data and handle all of the plotting for them.
    
    Plots:
        - Bar plot of sequence counts per cancer type.
        - Pie chart of sequences per cancer type in the whole dataset.
        - Lineplot of sequence length distribution per cancer type.
        - Heatmap of gene usage per cancer type.
        - Heatmap of Jaccard index between cancer types.
        - Heatmap of proportional amino acid presence per cancer type.

    Attributes:
        logger: Logging instance for the class.
        eda: ExploratoryAnalyzer instance containing exploratory data analysis data.
        unsupervised_vis: UnsupervisedVisualizer instance containing the PCA data.
        eda_save_path: Path in which we save the exploratory data analysis plots.
        unsupervised_save_path: Path in which we save the PCA plots.
    '''

    def __init__(self,eda: ExploratoryAnalyzer, unsupervised_vis: UnsupervisedVisualizer):

        self.logger = logging.getLogger("Plotter")
        self.eda = eda
        self.unsupervised_vis = unsupervised_vis
        self.eda_save_path = Path('Figures/EDA/')
        self.unsupervised_save_path = Path('Figures/Unsupervised/')
    
    def plot_sequence_count_bar(self):
        """ 
        Plots a bar chart with the number of sequences each disease type contains. 
        """

        self.logger.info("Plotting sequence count bar chart.")

        fig = plt.figure(figsize=(12,15))
        sns.barplot(x=list(self.eda.sequence_counts.keys()), y=list(self.eda.sequence_counts.values()))
        plt.xlabel("Cancer Types")
        plt.ylabel("Number of Sequences")
        plt.xticks(fontsize=10, rotation=45)
        plt.yticks(fontsize=10)
        plt.ticklabel_format(style='plain', axis='y')
        plt.title("Comparison of Total Sequence Count per Cancer Type")
        plt.savefig(self.eda_save_path / 'Sequence_Counts')
        self.logger.info("Sequence bar chart saved.")

    def plot_sequence_count_pie(self):
        """
        Plots a pie chart with the proportion of each disease type in the database. 
        """

        self.logger.info("Plotting sequence count pie chart.")

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
        plt.savefig(self.eda_save_path / 'Cancer_Type_Proportions')
        self.logger.info("Database proportion pie chart saved.")

    def plot_sequence_len_distribution(self):
        """ 
        Plots a line chart with the sequence length distribution of each repertoire. 
        """

        self.logger.info("Plotting sequence length distribution.")

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
        plt.savefig(self.eda_save_path / 'Amino_Acid_Distributions')
        self.logger.info("Sequence length distribution chart saved.")

    def plot_gene_usage(self):
        """ 
        Plots a heatmap of the gene permutations in each disease for V, D and J region.
        """

        self.logger.info("Plotting gene usage heatmap.")

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
        plt.savefig(self.eda_save_path / 'Gene Usage Heatmaps')

    def plot_jaccard_index(self):
        """
        Plots a heatmap of the Jaccard value between repertoires.
        """

        self.logger.info("Plotting Jaccard Index Heatmap")

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
        plt.savefig(self.eda_save_path / 'Jaccard_Heatmap')


    # Plot the amino acid distribution for each repertoire (heatmap).
    def plot_aa_distribtuion(self):
        """
        Plots the amino acid distribution for each repertoire in a heatmap.
        """

        self.logger.info("Plotting amino acid distribution heatmaps.")

        # Plot.
        keys = list(self.eda.aa_counts.keys())
        fig, ax = plt.subplots(3, 3, figsize=(30, 20), sharey='row')
        cbar_ax = fig.add_axes([.91, .3, .03, .4])

        for i in range(len(keys)):
            data = pd.DataFrame(self.eda.aa_counts[keys[i]]).fillna(0)
            sns.heatmap(data, vmax=0.4, cbar=i == 0, ax=ax.flat[i], cbar_ax=cbar_ax)
            ax.flat[i].set_title(f"{keys[i]} Amino Acid Distribution")
        fig.tight_layout(rect=[0, 0, .9, 1])
        plt.savefig(self.eda_save_path / 'Amino_Acid_Distributions')

    def make_pca_legend(self, data_dir:Path=Path('data/cleaned_files/')) -> dict[int, str]:
        """
        Returns the sorted list of diseases to match the sorted targets created in the Kmer Preprocessor.

        Args:
            data_dir: Path to directory where the function can access the target names. 

        Returns:
            Dictionary of {int: disease name} for each disease in the directory. 
        """
        return {value : disease for value, disease in enumerate(sorted(os.listdir(data_dir)))}
    
    def plot_pca(self, figure_title:str="PCA_scatter_plot"):
        """
        Plots a scatter graph of the pca data.
        """

        self.logger.info("Plotting PCA scatter plot.")

        # Data.
        categorical_labels = self.make_pca_legend()
        # Map the numerical targets to their categorical values for the legend. 
        pca_targets = self.unsupervised_vis.pca_targets.map(categorical_labels)

        # Plot.
        fig = plt.figure(figsize=(10, 10))
        sns.scatterplot(data=self.unsupervised_vis.pca_data, x=0, 
                        y=1, hue=pca_targets, palette='tab10', alpha=0.3)
        plt.title(figure_title)
        plt.xlabel("PCA 1")
        plt.ylabel("PCA 2")
        plt.savefig(self.unsupervised_save_path / f'{figure_title}.png')

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