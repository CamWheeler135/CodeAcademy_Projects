''' File that contains the custom errors. '''

class FileNotFoundError(Exception):
    ''' This error exists to raise errors if the user tries to collect data from an empty directory.'''

    def __str__(self) -> str:
        return f"No files found in the directory, ensure that the path name is correct or that files are loaded in the correct directory."
    

class DiseaseNotSupportedError(Exception):
    ''' This error exists to raise errors if users try to feed cancer type data into the pipeline that the model has not been trained on. '''

    def __str__(self, cancer_type) -> str:
        return f"{cancer_type} cancer is not in the supported list, of cancer types. Check filename or file before continuing."
    

class DataNotFoundError(Exception):
    """
    This error is raised when the data needed for EDA cannot be found.
    """

    def __str__(self) -> str:
        return f"The data needed for analysis are 'AASeq', 'Vregion', 'Dregion', 'Jregion', 'cloneFraction', please ensure they are all present."
