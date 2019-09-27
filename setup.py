from setuptools import setup, find_packages


setup(
    name='circmimi',
    version='0.3.0',
    url='https://github.com/TreesLab/CircMiMi',
    packages=find_packages(),
    install_requires=[
        'Click',
        'sqlalchemy',
        'pandas'
    ],
    entry_points={
        'console_scripts': [
            'circmimi_tools = circmimi.scripts.circmimi_tools:cli'
        ]
    }
)
