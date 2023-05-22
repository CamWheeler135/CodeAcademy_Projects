''' Base logger that all logs can be build off in each file. '''

# Module Imports.
import logging 
import sys
from rich.logging import RichHandler

# For logging to the log file. 
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%d-%m-%Y %H:%M',
                    filename='./logs.log',
                    filemode='w'
                    )

# For logging to the console
console_logger = RichHandler(markup=True)
console_handler = logging.StreamHandler(sys.stdout)
console_logger.setLevel(logging.INFO)
console_format = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
console_logger.setFormatter(console_format)
logging.getLogger('').addHandler(console_logger)
