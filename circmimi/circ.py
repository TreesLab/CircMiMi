import pandas as pd
import numpy as np
from collections import Counter, defaultdict
from circmimi.annotation import Annotator
from circmimi.ambiguous import AmbiguousChecker


class CircEvents:
    INPUT_COLUMNS = ('chr', 'pos1', 'pos2', 'strand', 'circ_id')
    DONOR_ACCEPTOR_COLUMNS = ('chr', 'donor', 'acceptor', 'strand', 'circ_id')

    def __init__(self, filename):
        self._filename = filename
        self.original_df = self._read_file(self._filename)
        self.df = self._get_donor_acceptor_df(self.original_df)

        self._summary_columns = {
            'summary': [],
            'filters': []
        }

    def _read_file(self, filename):
        df = pd.read_csv(filename, sep='\t', header=None, nrows=1)
        if df.shape[1] > 4:
            self.circ_ids_specified = True
        else:
            self.circ_ids_specified = False

        df = pd.read_csv(
            filename,
            sep='\t',
            names=self.INPUT_COLUMNS,
            dtype={
                'chr': 'category',
                'pos1': 'int',
                'pos2': 'int',
                'strand': 'category'
            }
        ).rename_axis('ev_id')

        return df

    @classmethod
    def _get_donor_acceptor(cls, s):
        chr_, pos1, pos2, strand, circ_id = s[list(cls.INPUT_COLUMNS)]
        pos1, pos2 = sorted([pos1, pos2])

        if strand == '+':
            donor_site = pos2
            acceptor_site = pos1
        elif strand == '-':
            donor_site = pos1
            acceptor_site = pos2

        res = [chr_, donor_site, acceptor_site, strand, circ_id]
        res = pd.Series(res, index=cls.DONOR_ACCEPTOR_COLUMNS)

        return res

    @classmethod
    def _get_donor_acceptor_df(cls, original_df):
        if original_df.empty:
            df = pd.DataFrame([], columns=cls.DONOR_ACCEPTOR_COLUMNS)
        else:
            df = original_df.apply(cls._get_donor_acceptor, axis=1)

        return df

    def expand_to_all_events(self, ev_df, fillna_value):
        expanded_df = self.original_df.reset_index()[['ev_id']].merge(
            ev_df,
            on='ev_id',
            how='left'
        ).fillna(fillna_value).set_index('ev_id')

        return expanded_df

    def submit_to_summary(self, summary_column, type_):
        self._summary_columns[type_].append(summary_column)

    @staticmethod
    def _get_host_genes(anno_df):

        def get_gene_symbol(anno_df):
            return anno_df['transcript'].apply(
                lambda t: t.gene.gene_symbol
            )

        host_gene_df = anno_df.assign(
            host_gene=get_gene_symbol
        )[[
            'ev_id',
            'host_gene'
        ]].drop_duplicates(
        ).groupby(
            'ev_id'
        ).agg(lambda genes: ','.join(sorted(genes))).reset_index()

        return host_gene_df

    @staticmethod
    def _get_circ_ids(host_genes):
        genes = host_genes['host_gene']
        gene_count = Counter(genes)
        idx_dict = defaultdict(int)

        circ_ids = []
        for gene in genes:
            if gene_count[gene] == 1:
                circ_ids.append(f'circ{gene}')
            else:
                idx_dict[gene] += 1
                circ_ids.append(f'circ{gene}_{idx_dict[gene]}')

        circ_ids_df = host_genes.assign(
            circ_id=pd.Series(circ_ids)
        )[['ev_id', 'circ_id']]

        return circ_ids_df

    def _update_circ_ids(self, circ_ids):
        all_ev_with_circ_id = self.expand_to_all_events(circ_ids, np.nan)
        self.original_df = self.original_df.assign(circ_id=all_ev_with_circ_id)
        self.df = self.df.assign(circ_id=all_ev_with_circ_id)

    def _append_host_genes(self, host_genes):
        all_ev_with_host_gene = self.expand_to_all_events(host_genes, np.nan)
        self.original_df = self.original_df.assign(host_gene=all_ev_with_host_gene)
        self.df = self.df.assign(host_gene=all_ev_with_host_gene)

    def check_annotation(self, anno_db_file):
        self._annotator = Annotator(anno_db_file)
        self.anno_df, anno_status = self.df.pipe(self._annotator.annotate)

        self._host_genes = self._get_host_genes(self.anno_df)
        self._append_host_genes(self._host_genes)

        if not self.circ_ids_specified:
            self._circ_ids = self._get_circ_ids(self._host_genes)
            self._update_circ_ids(self._circ_ids)

        self.submit_to_summary(anno_status, type_='filters')

    def check_ambiguous(self,
                        anno_db_file,
                        ref_file,
                        other_ref_file,
                        work_dir='.',
                        num_proc=1):

        self.checker = AmbiguousChecker(
            anno_db_file,
            ref_file,
            other_ref_file,
            work_dir=work_dir,
            num_proc=num_proc
        )
        checking_result = self.checker.check(self.df)

        self.submit_to_summary(checking_result, type_='filters')

    @staticmethod
    def _merge_columns(df, column_dfs):
        merged_df = df

        for col_df in column_dfs:
            merged_df = merged_df.merge(col_df, on='ev_id', how='left')

        return merged_df

    def get_filters_results(self, show_detail=False):
        filters_df = pd.DataFrame(self.df.index).pipe(
            self._merge_columns,
            self._summary_columns['filters']
        )

        if show_detail:
            filters_df = self.original_df.pipe(
                self._merge_columns,
                [filters_df]
            ).set_index('ev_id')

        return filters_df

    def get_summary(self):
        filters_df = self.get_filters_results()

        pass_column = filters_df.fillna('2').set_index(
            'ev_id'
        ).astype(
            int
        ).agg(
            'sum',
            axis=1
        ).eq(
            3
        ).replace(
            {
                True: 'yes',
                False: 'no'
            }
        ).reset_index().rename({0: 'pass'}, axis=1)

        summary_df = self.original_df.reset_index().pipe(
            self._merge_columns,
            self._summary_columns['summary'] + [pass_column, filters_df.fillna('NA')]
        ).set_index('ev_id')

        return summary_df

    @property
    def _passed_events(self):
        pass_df = self.get_summary()[lambda df: df['pass'] == 'yes'].reset_index()[['ev_id']]
        return pass_df

    @property
    def clear_df(self):
        return self.original_df.merge(self._passed_events, on='ev_id').set_index('ev_id')

    @property
    def clear_anno_df(self):
        return self.anno_df.merge(self._passed_events, on='ev_id', how='inner')

    @property
    def region_id(self):
        region_id_df = self.df.apply(
            lambda s: f"{s['chr']}:{s['donor']}|{s['acceptor']}({s['strand']})",
            axis=1
        ).rename('region_id')

        return region_id_df
