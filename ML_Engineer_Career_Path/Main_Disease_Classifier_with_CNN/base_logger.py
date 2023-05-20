''' Base logger that all logs can be build off in each file. '''

import logging 

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%d-%m-%Y %H:%M',
                    filename='./logs.log',
                    filemode='w'
                    )

