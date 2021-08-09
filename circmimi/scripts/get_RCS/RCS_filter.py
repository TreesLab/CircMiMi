#! /usr/bin/env python

import argparse
from collections import namedtuple


TITLES = [
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


class RCSReader:
    RCS_data = namedtuple('RCS_data', TITLES)

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


def RCS_filter(RCS_result,
               min_matches=80,
               min_alignment_len=50,
               min_bitscore=100):

    for data in RCS_result:
        if data.matches != '':
            if (float(data.matches) > min_matches) \
                and (int(data.alignment_len) > min_alignment_len) \
                    and (float(data.bit_score) > min_bitscore):

                yield data


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--min-matches', type=float, default=80)
    parser.add_argument('-l', '--min-aln-len', type=int, default=50)
    parser.add_argument('-b', '--min-bitscore', type=float, default=100)
    parser.add_argument('RCS_file', type=argparse.FileType('r'))

    return parser


def cli():
    parser = create_parser()
    args = parser.parse_args()

    print(*TITLES, sep='\t')

    reader = RCSReader(args.RCS_file)
    reader = RCS_filter(
        reader,
        min_matches=args.min_matches,
        min_alignment_len=args.min_aln_len,
        min_bitscore=args.min_bitscore
    )

    for data in reader:
        print(*data, sep='\t')


if __name__ == "__main__":
    cli()
