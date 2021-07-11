import logging


__version__ = '0.16.2'


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

ch = logging.StreamHandler()

formatter = logging.Formatter('%(asctime)s - %(message)s')
ch.setFormatter(formatter)

logger.addHandler(ch)
