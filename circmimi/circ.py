import pandas as pd
from circmimi.annotation import AnnotationUtils
from circmimi.ambiguous import AmbiguousChecker


class CircEvents:
    def __init__(self, filename):
        self._filename = filename
        self.original_df = self._read_file(self._filename)
        self.df = self._get_donor_acceptor_df(self.original_df)

        self._summary_columns = {
            'summary': [],
            'filters': [],
        }

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

    @classmethod
    def _get_donor_acceptor_df(cls, original_df):
        if original_df.empty:
            df = pd.DataFrame([], columns=['chr', 'donor', 'acceptor', 'strand'])
        else:
            df = original_df.apply(cls._get_donor_acceptor, axis=1)

        return df

    def check_annotation(self, anno_db):
        self._anno_utils = AnnotationUtils(anno_db)
        self.anno_df = self.df.pipe(self._anno_utils.get_annotation)

        no_common_transcript_df = self.anno_df.assign(
            no_common_transcript='0'
        )[[
            'ev_id',
            'no_common_transcript'
        ]].drop_duplicates(
        ).pipe(self._expand_to_all_events, fillna_value='1')

        self._submit_to_summary(no_common_transcript_df, type_='filters')

    def check_ambiguous(self,
                        ref_file,
                        other_ref_file,
                        work_dir='.',
                        num_proc=1):

        self.checker = AmbiguousChecker(
            ref_file,
            other_ref_file,
            work_dir=work_dir,
            num_proc=num_proc
        )
        self.checker.check(self.df)

        colinear_df = self.checker.colinear_result.assign(
            colinear='1'
        ).pipe(
            self._expand_to_all_events,
            fillna_value='0'
        )

        multiple_hits_df = self.checker.multiple_hits_result.assign(
            multiple_hits='1'
        ).pipe(
            self._expand_to_all_events,
            fillna_value='0'
        )

        self._submit_to_summary(colinear_df, type_='filters')
        self._submit_to_summary(multiple_hits_df, type_='filters')

    def _submit_to_summary(self, summary_column, type_):
        self._summary_columns[type_].append(summary_column)

    def _expand_to_all_events(self, ev_df, fillna_value):
        expanded_df = self.original_df.reset_index().rename(
            {
                'index': 'ev_id'
            },
            axis=1
        ).loc[:, ['ev_id']].merge(
            ev_df,
            on='ev_id',
            how='left'
        ).fillna(fillna_value)

        return expanded_df

    @staticmethod
    def _merge_columns(df, column_dfs):
        merged_df = df.reset_index().rename(
            {
                'index': 'ev_id'
            },
            axis=1
        )

        for col_df in column_dfs:
            merged_df = merged_df.merge(col_df, on='ev_id', how='left')

        merged_df = merged_df.drop('ev_id', axis=1)

        return merged_df

    def get_summary(self):
        summary_df = self.original_df.pipe(
            self._merge_columns,
            self._summary_columns['filters']
        )

        return summary_df

    @property
    def clear_df(self):
        return self.get_summary()[lambda df: df.iloc[:, 4:].agg('sum', axis=1) == 0].iloc[:, :4]

    @property
    def clear_anno_df(self):
        return self.anno_df.merge(
            self.clear_df.reset_index(),
            left_on='ev_id',
            right_on='index',
            how='inner'
        )
