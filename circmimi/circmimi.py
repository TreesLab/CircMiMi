import pandas as pd
from circmimi.annotation import AnnotationUtils
from circmimi.bed import BedUtils
from circmimi.seq import Seq
from circmimi.miranda import get_binding_sites, MirandaUtils


class Circmimi:
    @staticmethod
    def run(circ_file, anno_db_file, ref_file, mir_ref_file, num_proc=1):
        anno_utils = AnnotationUtils(anno_db_file)

        circ_events = CircEvents(circ_file)

        anno_df = circ_events.df.pipe(anno_utils.get_annotation)
        uniq_exons_df = anno_df.pipe(anno_utils.get_uniq_exons)
        bed_df = uniq_exons_df.pipe(BedUtils.to_bed)

        seq_df = bed_df.pipe(Seq.get_seq, ref_file=ref_file)
        seq_df['seq'] = seq_df.apply(Seq.extend_seq_for_circ_js, axis=1)

        miranda_df = seq_df.pipe(
            get_binding_sites,
            mir_ref_file=mir_ref_file,
            num_proc=num_proc
        ).pipe(
            MirandaUtils.append_exons_len,
            exons_len_df=uniq_exons_df[['exons_id', 'total_len']]
        ).pipe(
            MirandaUtils.remove_redundant_result
        ).pipe(
            MirandaUtils.append_cross_boundary
        ).pipe(
            MirandaUtils.append_ev_id,
            exons_ev_id_df=uniq_exons_df[['exons_id', 'ev_id']]
        ).pipe(
            MirandaUtils.append_merged_aln
        ).pipe(
            MirandaUtils.generate_aln_id
        ).drop('aln', axis=1)

        grouped_res_df = MirandaUtils.get_grouped_results(miranda_df)

        return grouped_res_df


class CircmimiResult:
    def __init__(self):
        pass


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
