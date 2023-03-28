import os
import gzip
import zipfile
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


def read_zip_file(zipped_file, member_file):
    with zipfile.ZipFile(zipped_file) as zip_files:
        with zip_files.open(member_file) as f_in:
            for line in f_in:
                yield line.decode()
