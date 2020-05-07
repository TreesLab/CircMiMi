import pandas as pd
from circmimi.annotation import AnnotationUtils


class CircEvents:
    def __init__(self, filename):
        self._filename = filename
        self.original_df = self._read_file(self._filename)
        self.df = self._get_donor_acceptor_df(self.original_df)

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

    @staticmethod
    def _append_status_column(df, status_df):
        status_df = df.reset_index().loc[:, ['index']].rename(
            {
                'index': 'ev_id'
            },
            axis=1
        ).merge(
            status_df,
            on='ev_id',
            how='left'
        ).drop(
            columns=['ev_id']
        ).fillna('1')

        df_with_status = pd.concat([df, status_df], axis=1)

        return df_with_status

    @property
    def status(self):
        no_common_transcript_df = self.anno_df.assign(
            no_common_transcript='0'
        )[[
            'ev_id',
            'no_common_transcript'
        ]].drop_duplicates()

        status_summary_df = self.original_df.pipe(
            self._append_status_column,
            status_df=no_common_transcript_df
        )

        return status_summary_df

    @property
    def clear_df(self):
        return self.status[lambda df: df.iloc[:, 4:].agg('sum', axis=1)==0]
