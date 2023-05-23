# Module Imports.
from pathlib import Path

# SRC Imports.
from src.data_cleaning import DataCleaner
from src.eda import ExploratoryAnalyzer

def main():
    cleaner = DataCleaner()
    cleaner.complete_clean()


if __name__ == "__main__":
    main()




