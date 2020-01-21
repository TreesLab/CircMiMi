import os
import gzip
from contextlib import contextmanager


@contextmanager
def cwd(path):
    origin_pwd = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(origin_pwd)


def open_file(filename):
    if filename.endswith('.gz'):
        opened_file = gzip.open(filename, 'rt')
    else:
        opened_file = open(filename)

    return opened_file
