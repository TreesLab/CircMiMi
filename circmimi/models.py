from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Index, Column, Integer, String, ForeignKey
from sqlalchemy.orm import relationship


Base = declarative_base()


class Chromosome(Base):
    __tablename__ = 'chromosome'

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return "Chromosome(name={})".format(self.name)


class Strand(Base):
    __tablename__ = 'strand'

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return "Strand(name={})".format(self.name)


class Biotype(Base):
    __tablename__ = 'biotype'

    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return 'Biotype(type="{}")'.format(self.name)


class DonorSite(Base):
    __tablename__ = 'donor_site'

    id = Column(Integer, primary_key=True)
    chr_id = Column(Integer, ForeignKey('chromosome.id'))
    junc_site = Column(Integer)
    strand_id = Column(Integer, ForeignKey('strand.id'))
    chr_ = relationship("Chromosome")
    strand = relationship("Strand")
    exons = relationship("Exon", back_populates='donor')

    __table_args__ = (Index('donor_site_index', "chr_id", "junc_site"), )

    def __init__(self, chr_id, junc_site, strand_id):
        self.chr_id = chr_id
        self.junc_site = junc_site
        self.strand_id = strand_id

    def __repr__(self):
        return "DonorSite({}, {}, {})".format(
            self.chr_.name, self.junc_site, self.strand.name
        )

    @property
    def region(self):
        return [self.chr_.name, self.junc_site, self.strand.name]


class AcceptorSite(Base):
    __tablename__ = 'acceptor_site'

    id = Column(Integer, primary_key=True)
    chr_id = Column(Integer, ForeignKey('chromosome.id'))
    junc_site = Column(Integer)
    strand_id = Column(Integer, ForeignKey('strand.id'))
    chr_ = relationship("Chromosome")
    strand = relationship("Strand")
    exons = relationship("Exon", back_populates='acceptor')

    __table_args__ = (Index('acceptor_site_index', "chr_id", "junc_site"), )

    def __init__(self, chr_id, junc_site, strand_id):
        self.chr_id = chr_id
        self.junc_site = junc_site
        self.strand_id = strand_id

    def __repr__(self):
        return "AcceptorSite({}, {}, {})".format(
            self.chr_.name, self.junc_site, self.strand.name
        )

    @property
    def region(self):
        return [self.chr_.name, self.junc_site, self.strand.name]


class Exon(Base):
    __tablename__ = 'exon'

    id = Column(Integer, primary_key=True)
    chr_id = Column(Integer, ForeignKey('chromosome.id'))
    start = Column(Integer)
    end = Column(Integer)
    strand_id = Column(Integer, ForeignKey('strand.id'))
    donor_id = Column(Integer, ForeignKey('donor_site.id'), index=True)
    acceptor_id = Column(Integer, ForeignKey('acceptor_site.id'), index=True)

    chr_ = relationship("Chromosome")
    strand = relationship("Strand")
    donor = relationship("DonorSite", back_populates='exons')
    acceptor = relationship("AcceptorSite", back_populates='exons')
    _transcript_exon = relationship('TranscriptExon', back_populates='exon')

    def __init__(self, chr_id, start, end, strand_id, donor_id, acceptor_id):
        self.chr_id = chr_id
        self.start = start
        self.end = end
        self.strand_id = strand_id
        self.donor_id = donor_id
        self.acceptor_id = acceptor_id

    def __repr__(self):
        return "Exon({}, {}, {}, {})".format(
            self.chr_.name, self.start, self.end, self.strand.name
        )

    def __len__(self):
        return self.end - self.start + 1

    @property
    def transcript(self):
        return list(map(lambda t_e: t_e.transcript, self._transcript_exon))

    @property
    def region(self):
        return [self.chr_.name, self.start, self.end, self.strand.name]


class Gene(Base):
    __tablename__ = 'gene'

    id = Column(Integer, primary_key=True)
    gene_id = Column(String)
    gene_symbol = Column(String)
    gene_type_id = Column(Integer, ForeignKey('biotype.id'))
    gene_type = relationship('Biotype')
    transcripts = relationship('Transcript', back_populates='gene')

    def __init__(self, gene_id, gene_symbol, gene_type_id):
        self.gene_id = gene_id
        self.gene_symbol = gene_symbol
        self.gene_type_id = gene_type_id

    def __repr__(self):
        return 'Gene(gene_id="{}", gene_symbol="{}", gene_type="{}")'.format(
            self.gene_id, self.gene_symbol, self.gene_type.name
        )


class Transcript(Base):
    __tablename__ = 'transcript'

    id = Column(Integer, primary_key=True)
    transcript_id = Column(String)
    gid = Column(Integer, ForeignKey('gene.id'))
    transcript_type_id = Column(Integer, ForeignKey('biotype.id'))
    transcript_type = relationship('Biotype')
    gene = relationship('Gene', back_populates='transcripts')
    _transcript_exon = relationship('TranscriptExon',
                                    back_populates='transcript')

    def __init__(self, transcript_id, gid, transcript_type_id):
        self.transcript_id = transcript_id
        self.gid = gid
        self.transcript_type_id = transcript_type_id

    def __repr__(self):
        return 'Transcript(transcript_id="{}", biotype="{}")'.format(
            self.transcript_id, self.transcript_type.name
        )

    def init_exons(self):
        self.exons = list(map(lambda t_e: t_e.exon,
                              sorted(self._transcript_exon,
                                     key=lambda t_e: t_e.exon_number)))

    @staticmethod
    def _get_len_of_exons(exons):
        return sum(map(len, exons))

    def get_inter_exons(self, from_, to_):
        from_, to_ = sorted([from_, to_])
        return self.exons[(from_ - 1):((to_ - 1) + 1)]

    def get_len_of_inter_exons(self, from_, to_):
        return self._get_len_of_exons(self.get_inter_exons(from_, to_))

    @property
    def total_len_of_exons(self):
        return self._get_len_of_exons(self.exons)

    @property
    def total_number_of_exons(self):
        return len(self.exons)

    def get_nth_exon(self, n):
        if (n >= 1) and (n <= self.total_number_of_exons):
            return self.exons[n - 1]
        else:
            return None

    def get_exon_number(self, exon):
        return self.exons.index(exon) + 1

    def get_len_of_intron(self, exon_num_1, up_or_down, empty_value='NA'):
        if up_or_down == 'up':
            exon_num_2 = exon_num_1 - 1
        elif up_or_down == 'down':
            exon_num_2 = exon_num_1 + 1

        exon_1 = self.get_nth_exon(exon_num_1)
        exon_2 = self.get_nth_exon(exon_num_2)
        if exon_1 and exon_2:
            junc_sites = sorted([exon_1.start, exon_1.end,
                                 exon_2.start, exon_2.end])
            intron_start = junc_sites[1] + 1
            intron_end = junc_sites[2] - 1
            return intron_end - intron_start + 1
        else:
            return empty_value

    @property
    def is_protein_coding(self):
        return self.transcript_type.name == 'protein_coding'

    def _get_exon_by_junc_site(self, donor_or_acceptor):
        if isinstance(donor_or_acceptor, DonorSite):
            junc_site_type = 'donor'
        elif isinstance(donor_or_acceptor, AcceptorSite):
            junc_site_type = 'acceptor'

        filtered_exons = list(filter(
            lambda exon: getattr(exon, junc_site_type) == donor_or_acceptor,
            self.exons
        ))

        if len(filtered_exons) == 1:
            return filtered_exons[0]
        elif len(filtered_exons) == 0:
            raise Exception("No exon founded!")
        else:
            raise Exception(
                "More than one exon was founded! Please check your database!"
            )

    def get_exon_by_donor(self, donor):
        return self._get_exon_by_junc_site(donor)

    def get_exon_by_acceptor(self, acceptor):
        return self._get_exon_by_junc_site(acceptor)


class TranscriptExon(Base):
    __tablename__ = 'transcript_exon'

    tid = Column(Integer, ForeignKey('transcript.id'),
                 primary_key=True, index=True)
    exon_number = Column(Integer, primary_key=True)
    eid = Column(Integer, ForeignKey('exon.id'), index=True)

    transcript = relationship('Transcript', back_populates='_transcript_exon')
    exon = relationship('Exon', back_populates='_transcript_exon')

    def __init__(self, tid, exon_number, eid):
        self.tid = tid
        self.exon_number = exon_number
        self.eid = eid

    def __repr__(self):
        return "{} {} {}".format(
            self.transcript.transcript_id, self.exon_number, self.exon
        )
