''' File that contains the code for cleaning data. '''


class DataCleaner():
    '''
    This class is for cleaning and standardizing the data for use in the ML pipeline. 
    - Collects the raw data from the 'raw_data' file.
    - Cleans the data, removing sequences based on length and read quality. 
    - Returns cleaned CSV files to 'cleaned_files'.
    '''
    
    def __init__(self):
        pass

    def collect_data(self):
        ''' Collect data from original files. '''
        pass

    def assert_size(self):
        ''' Remove sequences that exceed biological norms. '''
        pass

    def remove_bad_seqs(self):
        ''' Remove sequences that contain non-alphabetical characters. '''
        pass

    def name_files(self):
        ''' Gives files appropriate names according to disease. '''