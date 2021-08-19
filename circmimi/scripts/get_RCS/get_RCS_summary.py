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
    '#RCS across flanking sequences',
    '#RCS within the flanking sequence (the donor side)',
    '#RCS within the flanking sequence (the acceptor side)',
    '#RCS_across-#RCS_within>=1 (yes: 1; no: 0)'
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


def circRNA_id_reader(circ_file):
    if circ_file is not None:
        for line in circ_file:
            data = line.rstrip('\n').split('\t')
            chr_ = data[0]
            pos1, pos2 = sorted(map(int, data[1:3]))
            strand = data[3]

            if len(data) >= 5:
                circ_id = data[4]
            else:
                if strand == '+':
                    circ_id = f"{chr_}:{pos2}|{pos1}({strand})"
                elif strand == '-':
                    circ_id = f"{chr_}:{pos1}|{pos2}({strand})"

            yield circ_id


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--circ_file', type=argparse.FileType('r'),
                        help="4-columns TSV file: (chr, pos1, pos2, strand)")
    parser.add_argument('RCS_file', type=argparse.FileType('r'))

    return parser


def cli():
    parser = create_parser()
    args = parser.parse_args()

    print(*OUTPUT_TITLES, sep='\t')

    reader = RCSReader(args.RCS_file)
    group_reader = groupby(reader, key=lambda data: data.circRNA_id)

    circRNA_ids_reader = circRNA_id_reader(args.circ_file)

    for circRNA_id, RCS_gp in group_reader:

        for circ_id in circRNA_ids_reader:
            if circ_id != circRNA_id:
                print(circ_id, 0, 0, 0, sep='\t')
            else:
                break

        RCS_gp = list(RCS_gp)
        RCS_types = Counter(
            [(data.region1_type, data.region2_type) for data in RCS_gp]
        )

        across = RCS_types[('donor', 'acceptor')]
        within_donor = RCS_types[('donor', 'donor')]
        within_acceptor = RCS_types[('acceptor', 'acceptor')]
        across_minus_within = across - (within_donor + within_acceptor)
        across_minus_within_le_one = int(across_minus_within >= 1)

        print(
            circRNA_id,
            across,
            within_donor,
            within_acceptor,
            across_minus_within_le_one,
            sep='\t'
        )
    else:
        for circ_id in circRNA_ids_reader:
            print(circ_id, 0, 0, 0, 0, sep='\t')


if __name__ == "__main__":
    cli()
