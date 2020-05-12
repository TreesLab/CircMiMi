import pandas as pd
import numpy as np
import tempfile as tp
from functools import reduce
from circmimi.bed import BedUtils
from circmimi.seq import Seq
from circmimi.blat import Blat, PslFilters, PslUtils


class AmbiguousChecker:
    def __init__(self,
                 ref_file,
                 other_ref_file,
                 work_dir='.',
                 num_proc=1,
                 blat_bin='blat',
                 mp_blat_bin='mp_blat.py'):
        self.work_dir = work_dir
        self.num_proc = num_proc

        self.blat_bin = blat_bin
        self.mp_blat_bin = mp_blat_bin

        self.ref_file = ref_file
        self.other_ref_file = other_ref_file

    def check(self, circ_df):

        # 1. get flanking sequences
        flanking_regions = circ_df.apply(
            self._get_flanking_region,
            axis=1
        )
        flanking_regions_df = self._to_regions_df(flanking_regions)
        flanking_seq_df = flanking_regions_df.pipe(
            BedUtils.to_bed_df
        ).pipe(
            Seq.get_seq,
            ref_file=self.ref_file
        )

        fa_file = tp.NamedTemporaryFile(dir=self.work_dir)
        fa_file.write(Seq.to_fasta(flanking_seq_df).encode('utf-8'))
        fa_file.seek(0)

        # 2. blat
        blat = Blat(
            work_dir=self.work_dir,
            num_proc=self.num_proc,
            blat_bin=self.blat_bin,
            mp_blat_bin=self.mp_blat_bin
        )
        web_blat = Blat(
            work_dir=self.work_dir,
            num_proc=self.num_proc,
            blat_opts="-tileSize=9 -stepSize=9 -repMatch=32768",
            blat_bin=self.blat_bin,
            mp_blat_bin=self.mp_blat_bin
        )

        self._psl_rG_1 = blat(self.ref_file, fa_file.name)
        self._psl_rG_2 = web_blat(self.ref_file, fa_file.name)
        self._psl_rO_1 = blat(self.other_ref_file, fa_file.name)
        self._psl_rO_2 = web_blat(self.other_ref_file, fa_file.name)

        fa_file.close()

        # 3. evaluate
        # 3.1 colinear part
        psl_rG_1_colinear_df = self._psl_rG_1.df.pipe(PslFilters.colinear_filter)
        psl_rG_2_colinear_df = self._psl_rG_2.df.pipe(PslFilters.colinear_filter)
        psl_rO_1_colinear_df = self._psl_rO_1.df.pipe(PslFilters.colinear_filter)
        psl_rO_2_colinear_df = self._psl_rO_2.df.pipe(PslFilters.colinear_filter)

        psl_rG_1_colinear_ids = PslUtils.get_uniq_qname(psl_rG_1_colinear_df)
        psl_rG_2_colinear_ids = PslUtils.get_uniq_qname(psl_rG_2_colinear_df)
        psl_rO_1_colinear_ids = PslUtils.get_uniq_qname(psl_rO_1_colinear_df)
        psl_rO_2_colinear_ids = PslUtils.get_uniq_qname(psl_rO_2_colinear_df)

        self._all_colinear_ids = reduce(
            np.union1d,
            [
                psl_rG_1_colinear_ids,
                psl_rG_2_colinear_ids,
                psl_rO_1_colinear_ids,
                psl_rO_2_colinear_ids
            ]
        )

        # 3.2 multiple hits
        psl_rG_1_chimera_df = self._psl_rG_1.df.pipe(PslFilters.chimera_filter)
        psl_rG_2_chimera_df = self._psl_rG_2.df.pipe(PslFilters.chimera_filter)

        self._psl_rG_1_chimera_ids = PslUtils.get_uniq_qname(psl_rG_1_chimera_df)
        self._psl_rG_2_chimera_ids = PslUtils.get_uniq_qname(psl_rG_2_chimera_df)

        psl_rG_1_multiple_hits_df = self._psl_rG_1.df.pipe(
            PslUtils.remove_in_list,
            qnames=np.union1d(
                psl_rG_1_colinear_ids,
                self._psl_rG_1_chimera_ids
            )
        )

        psl_rG_2_multiple_hits_df = self._psl_rG_2.df.pipe(
            PslUtils.remove_in_list,
            qnames=np.union1d(
                psl_rG_2_colinear_ids,
                self._psl_rG_2_chimera_ids
            )
        )

        psl_rG_1_multiple_hits_ids = PslUtils.get_uniq_qname(psl_rG_1_multiple_hits_df)
        psl_rG_2_multiple_hits_ids = PslUtils.get_uniq_qname(psl_rG_2_multiple_hits_df)

        self._all_multiple_hits_ids = np.setdiff1d(
            np.union1d(
                psl_rG_1_multiple_hits_ids,
                psl_rG_2_multiple_hits_ids
            ),
            self._all_colinear_ids
        )

        self.colinear_result = pd.DataFrame(self._all_colinear_ids, columns=['ev_id'])
        self.multiple_hits_result = pd.DataFrame(self._all_multiple_hits_ids, columns=['ev_id'])

    @staticmethod
    def _get_flanking_region(circ_data):
        chr_ = circ_data.chr
        donor = circ_data.donor
        acceptor = circ_data.acceptor
        strand = circ_data.strand

        if strand == '+':
            donor_flanking = [chr_, donor - 99, donor, strand]
            acceptor_flanking = [chr_, acceptor, acceptor + 99, strand]
        elif strand == "-":
            donor_flanking = [chr_, donor, donor + 99, strand]
            acceptor_flanking = [chr_, acceptor - 99, acceptor, strand]

        return [donor_flanking, acceptor_flanking]

    @staticmethod
    def _to_regions_df(regions):
        df = pd.DataFrame(
            regions,
            columns=['regions']
        ).rename_axis('regions_id').reset_index()

        return df
