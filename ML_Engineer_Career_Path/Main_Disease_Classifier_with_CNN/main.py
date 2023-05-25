# Module Imports.
from pathlib import Path
import time

# SRC Imports.
from src.data_cleaning import DataCleaner
from src.eda import ExploratoryAnalyzer, Plotter

def main():
    # cleaner = DataCleaner()
    # cleaner.complete_clean()
    eda = ExploratoryAnalyzer()
    repertoires = eda.collect_files()
    eda.count_sequences_per_disease(repertoires)
    plotter = Plotter(eda)
    plotter.plot_sequence_count_pie()

if __name__ == "__main__":
    main()




