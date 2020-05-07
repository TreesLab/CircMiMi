import os.path
import pandas as pd
from circmimi.circ import CircEvents
from circmimi.bed import BedUtils
from circmimi.seq import Seq
from circmimi.miranda import get_binding_sites, MirandaUtils


class Circmimi:
    def __init__(self, work_dir='.'):
        self.work_dir = work_dir

        self.circ_events = None
        self.uniq_exons_df = None
        self.bed_df = None
        self.seq_df = None
        self.miranda_df = None
        self.grouped_res_df = None

    def run(self,
            circ_file,
            anno_db_file,
            ref_file,
            mir_ref_file,
            mir_target_file,
            num_proc=1):

        self.circ_events = CircEvents(circ_file)
        self.circ_events.check_annotation(anno_db_file)

        self.uniq_exons_df = self.circ_events.anno_df.pipe(self._get_uniq_exons)
        self.bed_df = self.uniq_exons_df.pipe(
            BedUtils.to_regions_df
        ).pipe(
            BedUtils.to_bed_df
        )

        self.seq_df = self.bed_df.pipe(Seq.get_extended_seq, ref_file=ref_file)

        self.miranda_df = self.seq_df.pipe(
            get_binding_sites,
            mir_ref_file=mir_ref_file,
            work_dir=self.work_dir,
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

        self.mir_target_db = get_mir_target_db(mir_target_file)

    def get_result_table(self):
        gene_symbol_df = self.circ_events.anno_df.assign(
            host_gene=lambda df: df['transcript'].apply(lambda t: t.gene.gene_symbol)
        ).loc[:, ['ev_id', 'host_gene']].drop_duplicates().reset_index(drop=True)

        res_df = self.circ_events.clear_df.reset_index().merge(
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
        ).merge(
            self.mir_target_db,
            on='mirna',
            how='left'
        )

        return res_df

    def save_circRNAs_status(self, out_file):
        if os.path.exists(out_file):
            status_df = pd.read_csv(out_file, sep='\t')
            status_df = status_df.merge(self.circ_events.status, on=['chr', 'pos1', 'pos2', 'strand'])
        else:
            status_df = self.circ_events.status

        status_df.to_csv(out_file, sep='\t', index=False)

    @staticmethod
    def _get_total_length(list_of_obj):
        return sum(map(len, list_of_obj))

    @classmethod
    def _get_uniq_exons(cls, anno_df):
        if anno_df.empty:
            uniq_exons_df = pd.DataFrame(
                [],
                columns=['exons', 'ev_id', 'total_len', 'exons_id']
            )
        else:
            uniq_exons_df = anno_df[['exons', 'ev_id']]\
                .drop_duplicates()\
                .reset_index(drop=True)

            uniq_exons_df['total_len'] = uniq_exons_df.apply(
                lambda s: cls._get_total_length(s['exons']),
                axis=1
            )

            uniq_exons_df['exons_id'] = uniq_exons_df.apply(
                lambda s: "exons_{}".format(s.name),
                axis=1
            )

        return uniq_exons_df


def get_mir_target_db(mir_tar_db_path):
    db = pd.read_csv(mir_tar_db_path, sep='\t', dtype='object')

    assert list(db.columns[:2]) == ['mirna', 'target_gene'], \
        ("The column names of the first two columns"
         " should be 'mirna' and 'target_gene'")

    return db
