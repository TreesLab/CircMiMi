import pandas as pd
import numpy as np
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from circmimi.models import Chromosome, Strand, DonorSite, AcceptorSite


class Annotation:
    def __init__(self, anno_db_file):
        self._anno_db_file = anno_db_file

        self.engine = create_engine('sqlite:///{}'.format(anno_db_file))
        self.Session = sessionmaker(bind=self.engine)
        self.session = self.Session()

        self.chr_dict = dict(
            self.session.query(Chromosome.name, Chromosome.id).all()
        )

        self.strand_dict = dict(
            self.session.query(Strand.name, Strand.id).all()
        )

    def _get_junc_site(self, JuncSiteType, chr_name, pos, strand=None):
        result = self.session\
            .query(JuncSiteType)\
            .filter_by(chr_id=self.chr_dict[chr_name], junc_site=pos)

        if strand is not None:
            result = result.filter_by(strand_id=self.strand_dict[strand])

        return result.first()

    def get_donor_site(self, chr_name, pos, strand=None):
        return self._get_junc_site(DonorSite, chr_name, pos, strand)

    def get_acceptor_site(self, chr_name, pos, strand=None):
        return self._get_junc_site(AcceptorSite, chr_name, pos, strand)

    @staticmethod
    def get_transcripts_of_exons(exons):
        transcripts = []
        for exon in exons:
            transcripts += exon.transcript

        return transcripts

    @classmethod
    def get_transcripts_of_js(cls, junc_site):
        if junc_site is None:
            return []

        exons = junc_site.exons
        transcripts = cls.get_transcripts_of_exons(exons)
        return transcripts

    @staticmethod
    def get_inter_exons_data(transcript, donor, acceptor):
        transcript.init_exons()

        exon_donor = transcript.get_exon_by_donor(donor)
        exon_acceptor = transcript.get_exon_by_acceptor(acceptor)

        exon_number_donor = transcript.get_exon_number(exon_donor)
        exon_number_acceptor = transcript.get_exon_number(exon_acceptor)

        inter_exons = tuple(
            transcript.get_inter_exons(exon_number_donor, exon_number_acceptor)
        )

        anno_data = [transcript,
                     exon_number_donor,
                     exon_number_acceptor,
                     inter_exons]

        return anno_data

    # @staticmethod
    # def get_exon_region(exon):
    #     return [exon.chr_.name, exon.start, exon.end, exon.strand.name]


class AnnotationUtils:
    def __init__(self, anno_db_file):
        self.db = Annotation(anno_db_file)

    def _get_anno_data_of_ev(self, s):
        chr_, donor_site, acceptor_site, strand = \
            s[['chr', 'donor', 'acceptor', 'strand']]

        ev_id = s.name

        donor = self.db.get_donor_site(chr_, donor_site, strand)
        acceptor = self.db.get_acceptor_site(chr_, acceptor_site, strand)

        donor_transcripts = self.db.get_transcripts_of_js(donor)
        acceptor_transcripts = self.db.get_transcripts_of_js(acceptor)

        common_transcripts = sorted(
            set(donor_transcripts) & set(acceptor_transcripts),
            key=lambda transcript: transcript.transcript_id
        )

        transcripts_data = [
            [ev_id] + self.db.get_inter_exons_data(transcript, donor, acceptor)
            for transcript in common_transcripts
        ]

        df_cols = [
            'ev_id',
            'transcript',
            'exon_number_d',
            'exon_number_a',
            'exons'
        ]

        transcripts_data_df = pd.DataFrame(transcripts_data, columns=df_cols)

        if not transcripts_data_df.empty:
            int_cols = [
                'ev_id',
                'exon_number_d',
                'exon_number_a'
            ]

            transcripts_data_df[int_cols] = \
                transcripts_data_df[int_cols].replace('NA', np.nan).astype('Int64')

        return transcripts_data_df

    def get_annotation(self, df):
        if df.empty:
            anno_df = pd.DataFrame(
                [],
                columns=[
                    'ev_id',
                    'transcript',
                    'exon_number_d',
                    'exon_number_a',
                    'exons'
                ]
            )
        else:
            raw_anno_dfs = df.apply(self._get_anno_data_of_ev, axis=1)
            anno_df = pd.concat(list(raw_anno_dfs)).reset_index(drop=True)

        return anno_df
