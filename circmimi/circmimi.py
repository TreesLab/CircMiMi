import pandas as pd
from circmimi.circ import CircEvents
from circmimi.bed import BedUtils
from circmimi.seq import Seq
from circmimi.miranda import get_binding_sites, MirandaUtils
from circmimi.rbp import PosMap, RBPBindingSites, RBPBindingSitesFilters


class Circmimi:
    def __init__(self,
                 anno_db_file,
                 ref_file,
                 mir_ref_file,
                 mir_target_file,
                 AGO_binding_file,
                 RBP_binding_file,
                 other_ref_file=None,
                 work_dir='.',
                 num_proc=1,
                 miranda_options=None):

        self.anno_db_file = anno_db_file
        self.ref_file = ref_file
        self.mir_ref_file = mir_ref_file
        self.mir_target_file = mir_target_file
        self.AGO_binding_file = AGO_binding_file
        self.RBP_binding_file = RBP_binding_file
        self.other_ref_file = other_ref_file
        self.work_dir = work_dir
        self.num_proc = num_proc
        self.miranda_options = miranda_options

        self.circ_events = None
        self.uniq_exons_df = None
        self.bed_df = None
        self.seq_df = None
        self.miranda_df = None
        self.grouped_res_df = None

        self.AGO_binding_sites = RBPBindingSites(self.AGO_binding_file)
        self.RBP_binding_sites = RBPBindingSites(self.RBP_binding_file)

    def run(self, circ_file):

        self.circ_events = CircEvents(circ_file)
        self.circ_events.check_annotation(self.anno_db_file)

        if self.other_ref_file is not None:
            self.circ_events.check_ambiguous(
                self.anno_db_file,
                self.ref_file,
                self.other_ref_file,
                work_dir=self.work_dir,
                num_proc=self.num_proc
            )

        self.uniq_exons_df = self.circ_events.clear_anno_df.pipe(
            self._get_uniq_exons
        )
        self.uniq_exons_regions_df = self.uniq_exons_df.pipe(
            BedUtils.to_regions_df
        )
        self.bed_df = self.uniq_exons_regions_df.pipe(
            BedUtils.to_bed_df
        )

        self.seq_df = self.bed_df.pipe(
            Seq.get_extended_seq,
            ref_file=self.ref_file
        )

        self.pos_map_db = {
            regions_id: PosMap(regions)
            for regions_id, regions in  self.uniq_exons_regions_df.values
        }

        # miRNAs part
        self.miranda_df = self.seq_df.pipe(
            get_binding_sites,
            mir_ref_file=self.mir_ref_file,
            work_dir=self.work_dir,
            num_proc=self.num_proc,
            miranda_options=self.miranda_options
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
            MirandaUtils.generate_uniq_id,
            column_name='aln'
        ).drop(
            'aln',
            axis=1
        ).pipe(
            MirandaUtils.get_genomic_position,
            pos_map_db=self.pos_map_db
        ).pipe(
            MirandaUtils.generate_uniq_id,
            column_name='genomic_regions'
        )

        self.miRNA_binding_sites_bed = self.miranda_df[[
            'genomic_regions_id',
            'genomic_regions'
        ]].drop_duplicates(
        ).reset_index(
            drop=True
        ).rename(
            {
                'genomic_regions_id': 'regions_id',
                'genomic_regions': 'regions'
            },
            axis=1
        ).pipe(
            BedUtils.to_bed_df,
            union=True
        )

        self.AGO_overlap_raw_data = self.AGO_binding_sites.overlap(
            self.miRNA_binding_sites_bed
        ).pipe(
            RBPBindingSites.append_joined_overlap
        )
        self.AGO_overlap = self.AGO_overlap_raw_data.pipe(
            RBPBindingSitesFilters.AGO_overlap_filter
        )

        self.AGO_overlap_count = self.AGO_overlap[[
            'name',
            'sample_id'
        ]].groupby(
            'name'
        ).count(
        ).reset_index(
        ).rename(
            {
                'name': 'genomic_regions_id',
                'sample_id': 'AGO_support'
            },
            axis=1
        ).astype(
            {
                'AGO_support': 'object'
            }
        )

        self.miranda_df = self.miranda_df.merge(
            self.AGO_overlap_count,
            on='genomic_regions_id',
            how='left'
        ).fillna(
            {
                'AGO_support': 0
            }
        ).assign(
            AGO_support_yn=lambda df: (df['AGO_support'] > 0).apply(int)
        )

        self.grouped_res_df = MirandaUtils.get_grouped_results(self.miranda_df)

        self.mir_target_db = get_mir_target_db(self.mir_target_file)

        self.gene_symbol_df = self.circ_events.clear_anno_df.assign(
            host_gene=lambda df: df['transcript'].apply(
                lambda t: t.gene.gene_symbol
            )
        )[['ev_id', 'host_gene']].drop_duplicates(
        ).sort_values(
            'host_gene'
        ).groupby(
            'ev_id'
        ).agg(','.join).reset_index()

        # final result table
        self.res_df = self.circ_events.clear_df.reset_index().merge(
            self.gene_symbol_df,
            on='ev_id',
            how='inner'
        ).merge(
            self.grouped_res_df,
            on='ev_id',
            how='inner'
        ).merge(
            self.mir_target_db,
            on='mirna',
            how='inner'
        ).sort_values(
            [
                'ev_id',
                'mirna',
                'target_gene'
            ]
        ).reset_index(drop=True)

        # submit summary
        circ_miRNA_count = self.res_df[['ev_id', 'mirna']].drop_duplicates().rename(
            {
                'mirna': '#circRNA_miRNA'
            },
            axis=1
        ).groupby(
            'ev_id'
        ).agg(
            'count'
        ).pipe(
            self.circ_events.expand_to_all_events,
            fillna_value=0
        ).astype('int')
        self.circ_events.submit_to_summary(circ_miRNA_count, type_='summary')

        circ_mRNA_count = self.res_df[['ev_id', 'target_gene']].drop_duplicates().rename(
            {
                'target_gene': '#circRNA_mRNA'
            },
            axis=1
        ).groupby(
            'ev_id'
        ).agg(
            'count'
        ).pipe(
            self.circ_events.expand_to_all_events,
            fillna_value=0
        ).astype('int')
        self.circ_events.submit_to_summary(circ_mRNA_count, type_='summary')

        total_count = self.res_df[['ev_id', 'mirna', 'target_gene']].drop_duplicates().assign(
            circRNA_miRNA_mRNA=1
        ).drop(
            ['mirna', 'target_gene'],
            axis=1
        ).groupby(
            'ev_id'
        ).agg(
            'count'
        ).pipe(
            self.circ_events.expand_to_all_events,
            fillna_value=0
        ).astype('int').rename(
            {
                'circRNA_miRNA_mRNA': '#circRNA_miRNA_mRNA'
            },
            axis=1
        )
        self.circ_events.submit_to_summary(total_count, type_='summary')

        self.res_df = self.res_df.drop('ev_id', axis=1)

    def save_result(self, out_file):
        self.res_df.to_csv(out_file, sep='\t', index=False)

    def save_circRNAs_summary(self, out_file):
        self.circ_events.get_summary().to_csv(out_file, sep='\t', index=False)

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
