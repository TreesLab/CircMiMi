import gzip
import re
from operator import itemgetter
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from circmimi.models import (Base, Chromosome, Strand, Biotype, Gene,
                             Transcript, Exon, TranscriptExon, DonorSite,
                             AcceptorSite)


class RegionAttr:
    attr_patter = re.compile(r'([^ ;]+) \"?([^;"]+)\"?;')

    def __init__(self, attr_string):
        self._attrs = dict(re.findall(self.attr_patter, attr_string))
        self._fix_attrs_from_ensembl()

    def _fix_attrs_from_ensembl(self):
        if "gene_biotype" in self._attrs:
            self._attrs["gene_type"] = self._attrs["gene_biotype"]

        if "transcript_biotype" in self._attrs:
            self._attrs["transcript_type"] = self._attrs["transcript_biotype"]

        if "gene_version" in self._attrs:
            self._attrs["gene_id"] = "{}.{}".format(
                self._attrs["gene_id"],
                self._attrs["gene_version"]
            )

        if "transcript_version" in self._attrs:
            self._attrs["transcript_id"] = "{}.{}".format(
                self._attrs["transcript_id"],
                self._attrs["transcript_version"]
            )

    def get(self, *attr_names):
        return [self._attrs.get(attr) for attr in attr_names]


class GetSite:
    get_start = itemgetter(0, 1, 3)
    get_end = itemgetter(0, 2, 3)

    @classmethod
    def get_donor_site(cls, exon):
        if exon[3] == 1:
            return cls.get_end(exon)
        elif exon[3] == 2:
            return cls.get_start(exon)

    @classmethod
    def get_acceptor_site(cls, exon):
        if exon[3] == 1:
            return cls.get_start(exon)
        elif exon[3] == 2:
            return cls.get_end(exon)


class TablesRawData:
    def __init__(self):
        pass

    def parse(self, anno_file):
        if anno_file.endswith('.gz'):
            opened_file = gzip.open(anno_file, 'rt')
        elif anno_file.endswith('.gtf'):
            opened_file = open(anno_file)
        else:
            raise Exception('File format not supported!')

        # Parse annotation gtf
        genes = []
        transcripts = []
        exons_tmp = []
        with opened_file as gz_in:
            for line in gz_in:
                if not line.startswith('#'):
                    data = line.rstrip('\n').split('\t')
                    region_type = data[2]
                    if region_type in ['gene', 'transcript', 'exon']:
                        attrs = RegionAttr(data[8])
                        if region_type == 'gene':
                            genes.append(
                                attrs.get(
                                    'gene_id',
                                    'gene_name',
                                    'gene_type'
                                )
                            )
                        elif region_type == 'transcript':
                            transcripts.append(
                                attrs.get(
                                    'transcript_id',
                                    'gene_id',
                                    'transcript_type'
                                )
                            )
                        elif region_type == 'exon':
                            exons_tmp.append(
                                list(itemgetter(0, 3, 4, 6)(data)) +
                                attrs.get(
                                    'transcript_id',
                                    'exon_number'
                                )
                            )

        chromosomes = sorted(set(map(itemgetter(0), exons_tmp)))
        strands = sorted(set(map(itemgetter(3), exons_tmp)))
        biotypes = sorted(set(map(itemgetter(2), genes + transcripts)))

        chromosomes_dict = self._get_index_map(chromosomes)
        strands_dict = self._get_index_map(strands)
        biotypes_dict = self._get_index_map(biotypes)
        genes_dict = self._get_index_map(map(itemgetter(0), genes))
        transcripts_dict = self._get_index_map(map(itemgetter(0), transcripts))

        # replace values by index number
        for g in genes:
            g[2] = biotypes_dict[g[2]]

        for t in transcripts:
            t[1] = genes_dict[t[1]]
            t[2] = biotypes_dict[t[2]]

        for e in exons_tmp:
            e[0] = chromosomes_dict[e[0]]
            e[1] = int(e[1])
            e[2] = int(e[2])
            e[3] = strands_dict[e[3]]
            e[4] = transcripts_dict[e[4]]
            e[5] = int(e[5])

        # get exons data
        exons = sorted(set(map(itemgetter(0, 1, 2, 3), exons_tmp)))
        exons_dict = self._get_index_map(exons)

        # get exon_transcript relation data
        exon_transcript = [[e[4], e[5], exons_dict[tuple(e[:4])]]
                           for e in exons_tmp]

        # get donor & acceptor
        donor_sites = sorted(set(map(GetSite.get_donor_site, exons)))
        acceptor_sites = sorted(set(map(GetSite.get_acceptor_site, exons)))

        donor_sites_dict = self._get_index_map(donor_sites)
        acceptor_sites_dict = self._get_index_map(acceptor_sites)

        # append the donor and acceptor indices to exons data
        exons_with_da = [list(e) +
                         [donor_sites_dict[GetSite.get_donor_site(e)],
                          acceptor_sites_dict[GetSite.get_acceptor_site(e)]]
                         for e in exons]

        self.chromosomes = chromosomes
        self.strands = strands
        self.biotypes = biotypes
        self.genes = genes
        self.transcripts = transcripts
        self.exons = exons_with_da
        self.exon_transcript = exon_transcript
        self.donor_sites = donor_sites
        self.acceptor_sites = acceptor_sites

    @staticmethod
    def _get_index_map(keys):
        return {k: i for i, k in enumerate(keys, start=1)}


def _write_data_to_db(session, data, DataModal):
    for row in data:
        if type(row) == str:
            session.add(DataModal(row))
        else:
            session.add(DataModal(*row))
    session.commit()


def generate(gtf_path, db_path):
    engine = create_engine('sqlite:///{}'.format(db_path))

    # create tables
    Base.metadata.create_all(bind=engine)

    Session = sessionmaker(bind=engine)
    session = Session()

    # parse raw data
    tables_raw_data = TablesRawData()
    tables_raw_data.parse(gtf_path)

    # write raw data to db
    _write_data_to_db(session, tables_raw_data.chromosomes, Chromosome)
    _write_data_to_db(session, tables_raw_data.strands, Strand)
    _write_data_to_db(session, tables_raw_data.biotypes, Biotype)
    _write_data_to_db(session, tables_raw_data.genes, Gene)
    _write_data_to_db(session, tables_raw_data.transcripts, Transcript)
    _write_data_to_db(session, tables_raw_data.exons, Exon)
    _write_data_to_db(session, tables_raw_data.exon_transcript, TranscriptExon)
    _write_data_to_db(session, tables_raw_data.donor_sites, DonorSite)
    _write_data_to_db(session, tables_raw_data.acceptor_sites, AcceptorSite)
