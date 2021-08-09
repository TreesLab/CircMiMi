#! /usr/bin/env python

import argparse
import subprocess as sp
import tempfile as tp
import io
import multiprocessing as mp
from collections import namedtuple
from collections.abc import Iterable
from operator import itemgetter


Bed6 = namedtuple('Bed6', ('chr_', 'start', 'end', 'name', 'score', 'strand'))


def get_fasta(bed_file, ref_file, bedtools_bin='bedtools'):
    cmd = [bedtools_bin, 'getfasta']
    cmd += ['-fi', ref_file]
    cmd += ['-bed', bed_file]
    cmd += ['-fo', '/dev/stdout']
    cmd += ['-s']
    # cmd += ['-name']

    res = sp.run(cmd, stdout=sp.PIPE, encoding='utf-8')

    return res.stdout


def blastn(query, subject, blastn_bin='blastn'):
    cmd = [blastn_bin]
    cmd += ['-task', 'blastn']
    cmd += ['-query', query]
    cmd += ['-subject', subject]
    cmd += ['-outfmt', '6']
    cmd += ['-strand', 'minus']
    # cmd += ['-word_size', '11']

    res = sp.run(cmd, stdout=sp.PIPE, encoding='utf-8')

    return res.stdout


def save_tmp(data, dir_='.'):
    tmp_file = tp.NamedTemporaryFile(dir=dir_)

    with open(tmp_file.name, 'w') as tmp_out:
        if isinstance(data, str):
            tmp_out.write(data)

        elif isinstance(data, Iterable):
            print(*data, sep='\t', file=tmp_out)

    return tmp_file


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


RCS_data = namedtuple('RCS_data', TITLES[1:13])


def read_blastn_result(blastn_result):
    for line in blastn_result.rstrip('\n').split('\n'):
        data = line.rstrip('\n').split('\t')
        data = RCS_data(*data)
        yield data


def is_disjoint(r1, r2):
    r1, r2 = sorted([r1, r2], key=itemgetter(0))

    if r1[1] < (r2[0] - 1):
        return True
    else:
        return False


def correct_results(blastn_result):
    all_alignments = []

    read_data = read_blastn_result(blastn_result)

    with io.StringIO() as out:
        for data in read_data:
            region1 = (int(data.alignment_start_1), int(data.alignment_end_1))
            region2 = (int(data.alignment_end_2), int(data.alignment_start_2))

            if is_disjoint(region1, region2):

                # check duplicates
                if (region2, region1) not in all_alignments:
                    all_alignments.append((region1, region2))

                    print(*data, sep='\t', file=out)

        return out.getvalue()


def read_circ_file(circ_file):
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

        if strand == '+':
            up_type = 'acceptor'
            down_type = 'donor'
        elif strand == '-':
            up_type = 'donor'
            down_type = 'acceptor'

        yield (chr_,
               pos1,
               pos2,
               strand,
               circ_id,
               up_type,
               down_type)


def get_RCS_output(circ_id, blastn_result, region1_type, region2_type):
    with io.StringIO() as tmp_out:
        if blastn_result:
            for res_line in blastn_result.rstrip('\n').split('\n'):
                print(circ_id, res_line, region1_type, region2_type, sep='\t', file=tmp_out)
        else:
            print(circ_id, *([''] * 12), region1_type, region2_type, sep='\t', file=tmp_out)

        return tmp_out.getvalue()


class RCS:
    def __init__(self, ref_file, dist, tmp_dir='.'):
        self.ref_file = ref_file
        self.dist = dist
        self.tmp_dir = tmp_dir

    def get_RCS(self, chr_, pos1, pos2, strand, circ_id, up_type, down_type):
        # upstream
        up_bed = Bed6(chr_, max(pos1 - self.dist, 0), pos1, circ_id, '.', strand)
        up_bed_file = save_tmp(up_bed, dir_=self.tmp_dir)

        up_fa = get_fasta(up_bed_file.name, self.ref_file)
        up_fa_file = save_tmp(up_fa, dir_=self.tmp_dir)

        # downstream
        down_bed = Bed6(chr_, pos2 - 1, pos2 - 1 + self.dist, circ_id, '.', strand)
        down_bed_file = save_tmp(down_bed, dir_=self.tmp_dir)

        down_fa = get_fasta(down_bed_file.name, self.ref_file)
        down_fa_file = save_tmp(down_fa, dir_=self.tmp_dir)

        # cross
        if strand == '+':
            blastn_result = blastn(down_fa_file.name, up_fa_file.name)
            RCS_result = get_RCS_output(circ_id, blastn_result, down_type, up_type)
        elif strand == '-':
            blastn_result = blastn(up_fa_file.name, down_fa_file.name)
            RCS_result = get_RCS_output(circ_id, blastn_result, up_type, down_type)

        # within upstream region
        blastn_result_up = blastn(up_fa_file.name, up_fa_file.name)
        if blastn_result_up:
            blastn_result_up = correct_results(blastn_result_up)
        RCS_result_up = get_RCS_output(circ_id, blastn_result_up, up_type, up_type)

        # within downstream region
        blastn_result_down = blastn(down_fa_file.name, down_fa_file.name)
        if blastn_result_down:
            blastn_result_down = correct_results(blastn_result_down)
        RCS_result_down = get_RCS_output(circ_id, blastn_result_down, down_type, down_type)

        return RCS_result, RCS_result_up, RCS_result_down


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dist', default=2000, type=int)
    parser.add_argument('ref_file')
    parser.add_argument('circ_file', type=argparse.FileType('r'), help="4-columns TSV file: (chr, pos1, pos2, strand)")
    parser.add_argument('-p', dest='num_proc', type=int, default=1)

    return parser


def cli():
    parser = create_parser()
    args = parser.parse_args()

    print(*TITLES, sep='\t')

    with tp.TemporaryDirectory(prefix='get_RCS.tmp.', dir='.') as tmp_dir:
        rcs = RCS(args.ref_file, args.dist, tmp_dir)

        with mp.Pool(processes=args.num_proc) as pool:
            RCS_results = [
                pool.apply_async(
                    rcs.get_RCS,
                    circ_data
                )
                for circ_data in read_circ_file(args.circ_file)
            ]

            for result in RCS_results:
                RCS_result, RCS_result_up, RCS_result_down = result.get()

                print(RCS_result, end='')
                print(RCS_result_up, end='')
                print(RCS_result_down, end='')


if __name__ == "__main__":
    cli()
