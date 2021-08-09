#! /usr/bin/env python

import argparse
from itertools import groupby
from collections import namedtuple, Counter


INPUT_TITLES = [
    'circRNA_id',
    'region1',
    'region2',
    'matches',
    'alignment_len',
    'mismatch',
    'gapopen',
    'alignment_start_1',
    'alignment_end_1',
    'alignment_start_2',
    'alignment_end_2',
    'E_value',
    'bit_score',
    'region1_type',
    'region2_type'
]


OUTPUT_TITLES = [
    'circRNA_id',
    '#RCS_across_flanking_introns',
    '#RCS_within_donor\'s_intron',
    '#RCS_within_acceptor\'s_intron'
]


class RCSReader:
    RCS_data = namedtuple('RCS_data', INPUT_TITLES)

    def __init__(self, RCS_result_file):
        self._RCS_result_file = RCS_result_file
        self._reader = self._RCS_data_reader()

    def _RCS_data_reader(self):
        self._title = self._RCS_result_file.readline()

        for line in self._RCS_result_file:
            data = line.rstrip('\n').split('\t')
            data = self.RCS_data(*data)
            yield data

    def __iter__(self):
        return self._reader


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('RCS_file', type=argparse.FileType('r'))

    return parser


def cli():
    parser = create_parser()
    args = parser.parse_args()

    print(*OUTPUT_TITLES, sep='\t')

    reader = RCSReader(args.RCS_file)
    group_reader = groupby(reader, key=lambda data: data.circRNA_id)

    for circRNA_id, RCS_gp in group_reader:
        RCS_gp = list(RCS_gp)
        RCS_types = Counter(
            [(data.region1_type, data.region2_type) for data in RCS_gp]
        )

        print(
            circRNA_id,
            RCS_types[('donor', 'acceptor')],
            RCS_types[('donor', 'donor')],
            RCS_types[('acceptor', 'acceptor')],
            sep='\t'
        )


if __name__ == "__main__":
    cli()
