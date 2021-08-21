import codecs
import os
import re
from setuptools import setup, find_packages


here = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    with codecs.open(os.path.join(here, *parts), 'r') as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


long_description = read('README.md')


setup(
    name='circmimi',
    version=find_version("circmimi", "__init__.py"),
    url='https://github.com/TreesLab/CircMiMi',
    packages=find_packages(),
    install_requires=[
        'click>=7.0',
        'sqlalchemy>=1.3.8',
        'numpy>=1.17.2',
        'pandas>=0.25.1',
        'openpyxl',
        'networkx>=2.4',
        'lxml>=4.5.0'
    ],
    entry_points={
        'console_scripts': [
            'circmimi_tools = circmimi.scripts.circmimi_tools:cli',
            'mp_blat.py = circmimi.scripts.mp_blat:cli',
            'get_RCS.py = circmimi.scripts.get_RCS.get_RCS:cli',
            'get_RCS_summary.py = circmimi.scripts.get_RCS.get_RCS_summary:cli',
            'RCS_filter.py = circmimi.scripts.get_RCS.RCS_filter:cli',
            'get_flanking_seq.py = circmimi.scripts.checkAA.get_flanking_seq:cli',
            'checkAA_reads.py = circmimi.scripts.checkAA.checkAA_reads:cli'
        ]
    },
    description=("A package for constructing CLIP-seq data-supported "
                 "\"circRNA - miRNA - mRNA\" interactions"),
    long_description=long_description,
    long_description_content_type="text/markdown"
)
