#! /usr/bin/env python

"""
A python script to check for reads if there exist colinear alignments
 or multiple alignments.

External tools:
  - blat
  - mp_blat.py

Usage:
  checkAA_reads.py -rg [Genome] -ro [OtherRefs] -p [NumProc] [Reads]
"""


import argparse
import tempfile as tp
import subprocess as sp
import io
import csv
import textwrap
import re
from itertools import groupby, chain
from functools import partial
from collections import namedtuple, deque


class Blat:
    def __init__(self, work_dir='.', blat_bin='blat', blat_opts="",
                 num_proc=1, mp_blat_bin='mp_blat.py'):
        self.work_dir = work_dir
        self.blat_bin = blat_bin
        self.blat_opts = blat_opts
        self.num_proc = num_proc
        self.mp_blat_bin = mp_blat_bin

    def _run(self, ref_file, fa_file):
        with tp.NamedTemporaryFile(dir=self.work_dir) as tmp_out:
            cmd = self._generate_cmd(ref_file, fa_file, tmp_out.name)
            sp.run(cmd)

            result = BlatPsl(tmp_out.read().decode('utf-8'))
            return result

    def _generate_cmd(self, ref_file, fa_file, out_file):
        cmd = [self.mp_blat_bin, ref_file, fa_file, out_file]
        cmd += ['-p', str(self.num_proc)]
        cmd += ['--tmp_path', self.work_dir]
        cmd += ['--blat_bin', self.blat_bin]
        cmd += ['--blat_opt', self.blat_opts]
        return cmd

    def __call__(self, ref_file, fa_file):
        return self._run(ref_file, fa_file)


class BlatPsl:
    PSL_TITLE = (
        'matches',
        'misMatches',
        'repMatches',
        'nCount',
        'qNumInsert',
        'qBaseInsert',
        'tNumInsert',
        'tBaseInsert',
        'strand',
        'qName',
        'qSize',
        'qStart',
        'qEnd',
        'tName',
        'tSize',
        'tStart',
        'tEnd',
        'blockCount',
        'blockSizes',
        'qStarts',
        'tStarts'
    )
    _str_cols = ['strand', 'qName', 'tName']
    _list_cols = ['blockSizes', 'qStarts', 'tStarts']
    _psl = namedtuple('Psl', PSL_TITLE)

    def __init__(self, psl_text):
        self.psl_text = psl_text
        self.psl_data = self._parse(psl_text)

    @classmethod
    def _parse(cls, psl_text):
        no_header_psl = cls._remove_header(psl_text)

        with io.StringIO() as tmp_out:
            print(*cls.PSL_TITLE, sep='\t', file=tmp_out)
            tmp_out.write(no_header_psl)
            tmp_out.seek(0)

            psl_data = [cls._psl(**cls._correct_column_type(data))
                        for data in csv.DictReader(tmp_out, delimiter='\t')]

            return psl_data

    @staticmethod
    def _remove_header(psl_text):
        if psl_text.startswith('psLayout'):
            no_header_psl = psl_text.split('\n', 5)[5]
        else:
            no_header_psl = psl_text

        return no_header_psl

    @classmethod
    def _correct_column_type(cls, data):
        corrected_data = {}
        for k, v in data.items():
            if k in cls._str_cols:
                new_v = v
            elif k in cls._list_cols:
                new_v = tuple(map(int, v.rstrip(',').split(',')))
            else:
                new_v = int(v)

            corrected_data[k] = new_v

        return corrected_data


class PslFilters:
    @staticmethod
    def _is_colinear(match, threshold=0.8):
        colinear_score = (match.matches + match.repMatches) / match.qSize
        return colinear_score > threshold

    @classmethod
    def colinear_filter(cls, psl_data, threshold=0.8):
        for match in psl_data:
            if cls._is_colinear(match, threshold=threshold):
                yield match

    @staticmethod
    def _is_chimera(psl_data, threshold):
        chimera = PslChimera()
        chimera.parse(psl_data)

        if not chimera.linear:
            if (threshold.check(chimera.five_end) is True) and \
               (threshold.check(chimera.three_end) is True):
                return True

        return False

    @classmethod
    def chimera_filter(cls, psl_data, min_matches=30, min_diff=3):
        threshold = _Threshold(min_matches=min_matches, min_diff=min_diff)
        _is_chimera = partial(cls._is_chimera, threshold=threshold)

        for _, gp in groupby(psl_data, key=lambda match: match.qName):
            gp = list(gp)
            if _is_chimera(gp):
                yield gp


class _Threshold:
    def __init__(self, min_matches, min_diff):
        self.min_matches = min_matches
        self.min_diff = min_diff

    def check(self, matches):
        num_matches = len(matches)
        if num_matches == 0:
            best_match_score = 0
            second_best_score = 0
        elif num_matches == 1:
            best_match_score = matches[0].matches
            second_best_score = 0
        elif num_matches == 2:
            best_match_score = matches[0].matches
            second_best_score = matches[1].matches

        if best_match_score < self.min_matches:
            return False

        if best_match_score - second_best_score < self.min_diff:
            return False

        return True


class PslChimera:
    def __init__(self):
        self._linear = False
        self._five_end = deque(maxlen=2)
        self._three_end = deque(maxlen=2)

    def parse(self, psl_data):
        for match in psl_data:
            if self._is_colinear(match):
                self._linear = True

            if self._is_five_end(match):
                self._add_best_match(self._five_end, match)

            if self._is_three_end(match):
                self._add_best_match(self._three_end, match)

    @property
    def linear(self):
        return self._linear

    @property
    def five_end(self):
        return self._five_end

    @property
    def three_end(self):
        return self._three_end

    @staticmethod
    def _is_colinear(match):
        return (match.qStart < 10) and (match.qSize - match.qEnd < 10)

    @staticmethod
    def _is_five_end(match):
        return match.qStart < 10

    @staticmethod
    def _is_three_end(match):
        return match.qSize - match.qEnd < 10

    @staticmethod
    def _is_better_match(match, other):
        return match.matches >= other.matches

    @classmethod
    def _add_best_match(cls, q, match):
        len_q = len(q)

        if len_q == 0:
            q.append(match)
        elif len_q == 1:
            if cls._is_better_match(match, q[0]):
                q.appendleft(match)
            else:
                q.append(match)
        else:
            if cls._is_better_match(match, q[0]):
                q.appendleft(match)
            elif cls._is_better_match(match, q[1]):
                q.appendleft(match)
                q.reverse()
            else:
                pass


class PslUtils:
    @staticmethod
    def get_uniq_qname(psl_data):
        uniq_qnames = sorted(set(map(lambda match: match.qName, psl_data)))
        return uniq_qnames


class AmbAlnChecker:
    def __init__(self,
                 ref_genome,
                 ref_others,
                 work_dir='.',
                 num_proc=1,
                 blat_bin='blat',
                 mp_blat_bin='mp_blat.py'):

        self.ref_genome = ref_genome
        self.ref_others = ref_others
        self.work_dir = work_dir
        self.num_proc = num_proc
        self.blat_bin = blat_bin
        self.mp_blat_bin = mp_blat_bin

        self.blat = Blat(
            work_dir=self.work_dir,
            num_proc=self.num_proc,
            blat_bin=self.blat_bin,
            mp_blat_bin=self.mp_blat_bin
        )
        self.web_blat = Blat(
            work_dir=self.work_dir,
            num_proc=self.num_proc,
            blat_opts="-tileSize=9 -stepSize=9 -repMatch=32768",
            blat_bin=self.blat_bin,
            mp_blat_bin=self.mp_blat_bin
        )

    @property
    def infos(self):
        infos = {
            'ref_genome': self.ref_genome,
            'ref_others': self.ref_others,
            'work_dir': self.work_dir,
            'num_proc': self.num_proc,
            'blat_bin': self.blat_bin,
            'mp_blat_bin': self.mp_blat_bin
        }
        return infos

    def check(self, reads_file):
        result = AmbAlnResult(reads_file=reads_file)
        result.set_infos(**self.infos)

        # run blat
        result._psl_rG_1 = self.blat(self.ref_genome, reads_file)
        result._psl_rG_2 = self.web_blat(self.ref_genome, reads_file)
        result._psl_rO_1 = self.blat(self.ref_others, reads_file)
        result._psl_rO_2 = self.web_blat(self.ref_others, reads_file)

        # evaluate
        # colinear part
        result._psl_rG_1_colinear = PslFilters.colinear_filter(result._psl_rG_1.psl_data)
        result._psl_rG_2_colinear = PslFilters.colinear_filter(result._psl_rG_2.psl_data)
        result._psl_rO_1_colinear = PslFilters.colinear_filter(result._psl_rO_1.psl_data)
        result._psl_rO_2_colinear = PslFilters.colinear_filter(result._psl_rO_2.psl_data)

        result._psl_rG_1_colinear_ids = PslUtils.get_uniq_qname(result._psl_rG_1_colinear)
        result._psl_rG_2_colinear_ids = PslUtils.get_uniq_qname(result._psl_rG_2_colinear)
        result._psl_rO_1_colinear_ids = PslUtils.get_uniq_qname(result._psl_rO_1_colinear)
        result._psl_rO_2_colinear_ids = PslUtils.get_uniq_qname(result._psl_rO_2_colinear)

        result.colinear_ids = sorted(set(chain(
            result._psl_rG_1_colinear_ids,
            result._psl_rG_2_colinear_ids,
            result._psl_rO_1_colinear_ids,
            result._psl_rO_2_colinear_ids
        )))

        # multiple hits part
        result._psl_rG_1_chimera = PslFilters.chimera_filter(result._psl_rG_1.psl_data)
        result._psl_rG_2_chimera = PslFilters.chimera_filter(result._psl_rG_2.psl_data)
        result._psl_rG_1_chimera_ids = PslUtils.get_uniq_qname(chain.from_iterable(result._psl_rG_1_chimera))
        result._psl_rG_2_chimera_ids = PslUtils.get_uniq_qname(chain.from_iterable(result._psl_rG_2_chimera))

        result._psl_rG_1_multiple_hits = list(filter(
            lambda m: m.qName not in result._psl_rG_1_colinear_ids + result._psl_rG_1_chimera_ids,
            result._psl_rG_1.psl_data
        ))
        result._psl_rG_2_multiple_hits = list(filter(
            lambda m: m.qName not in result._psl_rG_2_colinear_ids + result._psl_rG_2_chimera_ids,
            result._psl_rG_2.psl_data
        ))
        result._psl_rG_1_multiple_hits_ids = PslUtils.get_uniq_qname(result._psl_rG_1_multiple_hits)
        result._psl_rG_2_multiple_hits_ids = PslUtils.get_uniq_qname(result._psl_rG_2_multiple_hits)

        result.multiple_hits_ids = sorted(
            set(chain(
                result._psl_rG_1_multiple_hits_ids,
                result._psl_rG_2_multiple_hits_ids
            )).difference(
                set(result.colinear_ids)
            )
        )

        result.update_result()

        return result


class AmbAlnResult:
    _id_in_fasta = re.compile(r'^>([^ ]+).*\n$')

    def __init__(self, reads_file):
        self.reads_file = reads_file
        self.infos = {}

        self._init_result()

    def _init_result(self):
        self.read_ids = self._get_read_ids(self.reads_file)
        self.result = {read_id: [0, 0, 0] for read_id in self.read_ids}

    @classmethod
    def _get_read_ids(cls, reads_file):
        read_ids = []
        with open(reads_file) as f_in:
            for line in f_in:
                if line.startswith('>'):
                    read_ids.append(re.match(cls._id_in_fasta, line).group(1))
        return read_ids

    def set_infos(self, **kwargs):
        self.infos.update(kwargs)

    def update_result(self):
        for id_ in self.colinear_ids:
            self.result[id_][0] = 1

        for id_ in self.multiple_hits_ids:
            self.result[id_][1] = 1

        for id_, checking_result in self.result.items():
            if (checking_result[0] == 1) or (checking_result[1] == 1):
                self.result[id_][2] = 1

    def save_result(self, out_file):
        with open(out_file, 'w') as out:
            csv_writer = csv.writer(out, delimiter='\t')

            csv_writer.writerow(
                [
                    'circRNA_id',
                    'with an alternative co-linear explanation',
                    'with multiple_hits',
                    'alignment ambiguity (with an alternative co-linear explanation or multiple hits)'
                ]
            )

            for read_id, check_items in self.result.items():
                csv_writer.writerow([read_id] + check_items)


def create_parser():

    parser = argparse.ArgumentParser(

        description=textwrap.dedent("""
            A python script to check for reads if there exist colinear alignments
             or multiple alignments.
            """),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-rg', dest='ref_genome', required=True)
    parser.add_argument('-ro', dest='ref_others', required=True)
    parser.add_argument('in_file', help='input reads file.')
    parser.add_argument('out_file')
    parser.add_argument('-p', dest='num_proc', type=int, default=1)

    return parser


def cli():
    parser = create_parser()
    args = parser.parse_args()

    checker = AmbAlnChecker(
        ref_genome=args.ref_genome,
        ref_others=args.ref_others,
        num_proc=args.num_proc
    )

    result = checker.check(args.in_file)
    result.save_result(args.out_file)


if __name__ == "__main__":
    cli()
