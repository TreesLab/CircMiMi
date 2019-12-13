import io
import pandas as pd
import subprocess as sp
import tempfile as tp
from functools import partial
from collections import deque


class Blat:
    def __init__(self, work_dir='.', blat_bin='blat', blat_opts="", num_proc=1, mp_blat_bin='mp_blat.py'):
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

    def __init__(self, psl_text):
        self.psl_text = psl_text
        self.df = self.parse(psl_text)

    @classmethod
    def parse(cls, psl_text):
        no_header_psl = cls.remove_header(psl_text)

        with io.StringIO() as tmp_out:
            tmp_out.write(no_header_psl)
            tmp_out.seek(0)

            psl_df = pd.read_csv(
                tmp_out,
                sep='\t',
                header=None,
                names=cls.PSL_TITLE
            )

            return psl_df

    @staticmethod
    def remove_header(psl_text):
        if psl_text.startswith('psLayout'):
            no_header_psl = psl_text.split('\n', 5)[5]
        else:
            no_header_psl = psl_text

        return no_header_psl


class PslFilters:
    @staticmethod
    def _is_colinear(df, threshold=0.8):
        return ((df.matches + df.repMatches) / df.qSize) > threshold

    @classmethod
    def colinear_filter(cls, df, threshold=0.8):
        _is_colinear = partial(cls._is_colinear, threshold=threshold)
        return df[_is_colinear]

    @staticmethod
    def _is_chimera(subdf, threshold):
        chimera = PslChimera()
        chimera.parse(subdf)

        if not chimera.linear:
            if (threshold.check(chimera.five_end) is True) and \
               (threshold.check(chimera.three_end) is True):
                return True

        return False

    @classmethod
    def chimera_filter(cls, df, min_matches=30, min_diff=3):
        threshold = _Threshold(min_matches=min_matches, min_diff=min_diff)
        _is_chimera = partial(cls._is_chimera, threshold=threshold)
        chimera_df = df.groupby('qName').filter(_is_chimera)
        return chimera_df


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

    def parse(self, df):
        for idx, match in df.iterrows():
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
    def get_uniq_qname(df):
        return df.qName.unique()

    @staticmethod
    def remove_in_list(df, qnames):
        return df[~df.qName.isin(qnames)]
