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
        'Click',
        'sqlalchemy',
        'pandas',
        'xlrd'
    ],
    entry_points={
        'console_scripts': [
            'circmimi_tools = circmimi.scripts.circmimi_tools:cli'
        ]
    },
    description=("A toolset for investigating the interactions between"
                 " circRNA - miRNA - mRNA."),
    long_description=long_description
)
