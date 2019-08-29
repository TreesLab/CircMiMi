import pandas as pd
from circmimi.annotation import AnnotationUtils
from circmimi.bed import BedUtils
from circmimi.seq import Seq
from circmimi.miranda import get_binding_sites, MirandaUtils


class Circmimi:
    def __init__(self):
        self.circ_events = None
        self.anno_df = None
        self.uniq_exons_df = None
        self.bed_df = None
        self.seq_df = None
        self.miranda_df = None
        self.grouped_res_df = None

    def run(self, circ_file, anno_db_file, ref_file, mir_ref_file, num_proc=1):
        self.anno_utils = AnnotationUtils(anno_db_file)

        self.circ_events = CircEvents(circ_file)
        self.anno_df = self.circ_events.df.pipe(self.anno_utils.get_annotation)
        self.uniq_exons_df = self.anno_df.pipe(self.anno_utils.get_uniq_exons)
        self.bed_df = self.uniq_exons_df.pipe(BedUtils.to_bed)

        self.seq_df = self.bed_df.pipe(Seq.get_seq, ref_file=ref_file)
        self.seq_df['seq'] = self.seq_df.apply(
            Seq.extend_seq_for_circ_js,
            axis=1
        )

        self.miranda_df = self.seq_df.pipe(
            get_binding_sites,
            mir_ref_file=mir_ref_file,
            num_proc=num_proc
        ).pipe(
            MirandaUtils.append_exons_len,
            exons_len_df=self.uniq_exons_df[['exons_id', 'total_len']]
        ).pipe(
            MirandaUtils.remove_redundant_result
        ).pipe(
            MirandaUtils.append_cross_boundary
        ).pipe(
            MirandaUtils.append_ev_id,
            exons_ev_id_df=self.uniq_exons_df[['exons_id', 'ev_id']]
        ).pipe(
            MirandaUtils.append_merged_aln
        ).pipe(
            MirandaUtils.generate_aln_id
        ).drop('aln', axis=1)

        self.grouped_res_df = MirandaUtils.get_grouped_results(self.miranda_df)

    def get_result_table(self):
        gene_symbol_df = self.anno_df.apply(
            lambda s: pd.Series(
                [
                    s['ev_id'],
                    s['transcript'].gene.gene_symbol
                ],
                index=['ev_id', 'gene_symbol']
            ),
            axis=1
        ).drop_duplicates(
        ).reset_index(
            drop=True
        )

        res_df = self.circ_events.original_df.reset_index().merge(
            gene_symbol_df,
            left_on='index',
            right_on='ev_id',
            how='left'
        ).merge(
            self.grouped_res_df,
            on='ev_id',
            how='left'
        ).drop(
            ['index', 'ev_id'],
            axis=1
        )

        return res_df


class CircEvents:
    def __init__(self, filename):
        self._filename = filename
        self.original_df = self._read_file(self._filename)
        self.df = self.original_df.apply(self._get_donor_acceptor, axis=1)

    @staticmethod
    def _read_file(filename):
        df = pd.read_csv(
            filename,
            sep='\t',
            names=['chr', 'pos1', 'pos2', 'strand'],
            dtype={
                'chr': 'category',
                'pos1': 'int',
                'pos2': 'int',
                'strand': 'category'
            }
        )
        return df

    @staticmethod
    def _get_donor_acceptor(s):
        chr_, pos1, pos2, strand = s[['chr', 'pos1', 'pos2', 'strand']]
        pos1, pos2 = sorted([pos1, pos2])

        if strand == '+':
            donor_site = pos2
            acceptor_site = pos1
        elif strand == '-':
            donor_site = pos1
            acceptor_site = pos2

        res = [chr_, donor_site, acceptor_site, strand]
        res = pd.Series(res, index=['chr', 'donor', 'acceptor', 'strand'])

        return res
