import pandas as pd
import numpy as np
import tempfile as tp
from functools import reduce
from circmimi.bed import BedUtils
from circmimi.seq import Seq
from circmimi.blat import Blat, PslFilters, PslUtils
from circmimi.annotation import Annotation


class AmbiguousChecker:
    _CHECK_LIST = [
        'ambiguity_with_an_colinear_explanation',
        'ambiguity_with_multiple_hits'
    ]

    def __init__(self,
                 anno_db_file,
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

        self.anno_db = Annotation(anno_db_file)
        self.ref_file = ref_file
        self.other_ref_file = other_ref_file

    def check(self, circ_df):
        self._init_status(circ_df)

        # 1. get flanking sequences
        flanking_regions = circ_df.apply(
            self._get_flanking_region,
            axis=1
        )
        flanking_regions_df = self._to_regions_df(flanking_regions).dropna()
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

        for id_ in self._all_colinear_ids:
            self._report_status(id_, self._CHECK_LIST[0])

        for id_ in self._all_multiple_hits_ids:
            self._report_status(id_, self._CHECK_LIST[1])

        return self._checking_result

    def _get_flanking_region(self, circ_data):
        chr_ = circ_data.chr
        donor_site = circ_data.donor
        acceptor_site = circ_data.acceptor
        strand = circ_data.strand

        ev_id = circ_data.name

        donor = self.anno_db\
            .get_nearest_donor_site(chr_, donor_site, strand, dist=5)
        acceptor = self.anno_db\
            .get_nearest_acceptor_site(chr_, acceptor_site, strand, dist=5)

        if (donor is None) or (acceptor is None):
            self._report_status(ev_id, self._CHECK_LIST[0], np.nan)
            self._report_status(ev_id, self._CHECK_LIST[1], np.nan)
            return np.nan

        donor_site = donor.junc_site
        donor_exon = max(donor.exons, key=lambda exon: len(exon))
        donor_acceptor = donor_exon.acceptor.junc_site

        acceptor_site = acceptor.junc_site
        acceptor_exon = max(acceptor.exons, key=lambda exon: len(exon))
        acceptor_donor = acceptor_exon.donor.junc_site

        if strand == '+':

            donor_flanking = [
                chr_,
                max(donor_site - 99, donor_acceptor),
                donor_site,
                strand
            ]

            acceptor_flanking = [
                chr_,
                acceptor_site,
                min(acceptor_site + 99, acceptor_donor),
                strand
            ]

        elif strand == "-":
            donor_flanking = [
                chr_,
                donor_site,
                min(donor_site + 99, donor_acceptor),
                strand
            ]

            acceptor_flanking = [
                chr_,
                max(acceptor_site - 99, acceptor_donor),
                acceptor_site,
                strand
            ]

        return [donor_flanking, acceptor_flanking]

    @staticmethod
    def _to_regions_df(regions):
        df = pd.DataFrame(
            regions,
            columns=['regions']
        ).rename_axis('regions_id').reset_index()

        return df

    def _init_status(self, circ_df):
        self._checking_result = pd.DataFrame(
            [],
            columns=self._CHECK_LIST
        ).rename_axis('ev_id')

        for id_ in circ_df.index:
            self._report_status(id_)

    def _report_status(self, ev_id, status=None, value='1', init_value='0'):
        if ev_id in self._checking_result.index:
            if status is not None:
                self._checking_result.loc[ev_id, status] = value
        else:
            ev_status = pd.Series(
                [init_value] * len(self._CHECK_LIST),
                index=self._CHECK_LIST,
                name=ev_id
            )

            if status is not None:
                ev_status[status] = value

            self._checking_result = self._checking_result.append(ev_status)
