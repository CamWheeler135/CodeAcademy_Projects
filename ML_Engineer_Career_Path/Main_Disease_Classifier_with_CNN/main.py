# Module Imports.
from pathlib import Path
import time

# SRC Imports.
from src.data_cleaning import DataCleaner, DataHandler
from src.eda import ExploratoryAnalyzer, UnsupervisedVisualizer, Plotter
from src.preprocessing import UnsupervisedPreprocessor

def main():
    ## Data Cleaning. 
    dc_dh = DataHandler(collection_path=Path('data/raw_files/'), save_path=Path('data/cleaned_files/'))
    data_cleaner = DataCleaner(dc_dh)
    data_cleaner.complete_clean()

    ## Exploratory Data Analysis.
    eda_dh = DataHandler(collection_path=Path('data/cleaned_files'))
    eda = ExploratoryAnalyzer(eda_dh)
    eda.complete_eda()


    ## Unsupervised Preprocessing.
    usp_dh = DataHandler(
        collection_path=Path('data/cleaned_files'), 
        save_path=Path('data/processed_files/kmers/kmers.csv'))
    usp = UnsupervisedPreprocessor(usp_dh)
    usp.complete_preprocess()

    ## PCA Analysis
    pca_dh = DataHandler(collection_path=Path('data/processed_files/kmers'),
                         save_path=Path('data/processed_files/pca/pca_df.csv'))
    pca = UnsupervisedVisualizer(pca_dh)
    pca.complete_pca_analysis()

    ## Plotting
    plotter = Plotter(eda, pca)
    plotter.complete_plots()



if __name__ == "__main__":
    main()



