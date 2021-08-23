from collections import namedtuple


def read_table_with_titles(filename, sep='\t', titles=None):
    """A generator function for reading data.

    Args:
      filename: the path to the data file

    Yield:
      Data: a namedtuple with titles from the data
    """

    with open(filename) as f_in:

        if titles is None:
            titles = f_in.readline().rstrip('\n').split(sep)

        Data = namedtuple('Data', titles)

        for line in f_in:
            data = Data(*line.rstrip('\n').split(sep))
            yield data


def read_table(filename, sep='\t', header=True):
    with open(filename) as f_in:
        if header:
            _ = f_in.readline()

        for line in f_in:
            data = line.rstrip('\n').split(sep)
            yield data
