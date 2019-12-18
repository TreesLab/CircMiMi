import re
import os
import subprocess as sp
import pandas as pd
import csv
import shutil
from contextlib import contextmanager
from circmimi.reference import gendb
from circmimi.reference.species import species_list
from circmimi.seq import parse_fasta
from circmimi.reference import resource as rs


class RefFile:
    def __init__(self, src):
        self.src = src
        self.filename = None

    def generate(self):
        return self.filename


class AnnoRef(RefFile):
    def generate(self):
        self.filename = re.sub(r'\.gtf(?:\.gz)?$', '.db', self.src)
        gendb.generate(self.src, self.filename)
        return self.filename


class MirRef(RefFile):
    def generate(self, species):
        species = species.key

        unzipped_file = self.src.filename
        self.filename = re.sub(r'\.fa$', '.{}.fa'.format(species), unzipped_file)

        with open(unzipped_file) as fa_in:
            fa_txt = fa_in.read()

        with open(self.filename, 'w') as out:
            for fa_id, fa_seq in parse_fasta(fa_txt):
                m = re.search(r'^{}-'.format(species), fa_id)
                if m:
                    print(">{}\n{}".format(fa_id, fa_seq), file=out)

        return self.filename


class MiRTarBaseRef(RefFile):
    def generate(self, species):
        df = pd.read_excel(self.src)

        formatted_data = df[
            df['Species (miRNA)'] == species.fullname
        ][[
            'miRNA',
            'Target Gene',
            'References (PMID)'
        ]].groupby([
            'miRNA',
            'Target Gene'
        ]).agg(
            'count'
        ).reset_index(
        ).rename(
            {
                'miRNA': 'mirna',
                'Target Gene': 'target_gene',
                'References (PMID)': 'ref_count'
            },
            axis=1
        )

        self.filename = "miRTarBase.{}.tsv".format(species.key)

        formatted_data.to_csv(
            self.filename,
            sep='\t',
            index=False,
            quoting=csv.QUOTE_NONE
        )

        return self.filename


class EnsemblTranscriptsFile(RefFile):
    def generate(self, biotype):
        unzipped_file = self.src.filename
        self.filename = re.sub(r'\.fa$', '.{}.fa'.format(biotype), unzipped_file)

        with open(unzipped_file) as fa_in:
            fa_txt = fa_in.read()

        with open(self.filename, 'w') as out:
            for fa_id, fa_seq in parse_fasta(fa_txt):
                m = re.search(r'{}'.format(biotype), fa_id)
                if m:
                    print(">{}\n{}".format(fa_id, fa_seq), file=out)

        return self.filename


class RepChrM(RefFile):
    def generate(self):
        unzipped_file = self.src.filename
        file_path = os.path.dirname(unzipped_file)
        self.filename = os.path.join(file_path, "RepChrM.fa")

        with open(unzipped_file) as fa_in:
            fa_txt = fa_in.read()

        with open(self.filename, 'w') as out:
            for fa_id, fa_seq in parse_fasta(fa_txt):
                m = re.search(r'^(?:chr)?MT?', fa_id)
                if m:
                    print(">{}\n{}\n{}".format(fa_id, fa_seq, fa_seq),
                          file=out)

        return self.filename


class OtherTranscripts:
    def __init__(self, src_type, pc_src, lncRNA_src, repChrM_src):
        self.src_type = src_type
        self.pc_src = pc_src
        self.lncRNA_src = lncRNA_src
        self.repChrM_src = repChrM_src
        self.filename = 'others.fa'

    def generate(self):
        if self.src_type == 'ensembl':
            pc_ref = EnsemblTranscriptsFile(self.pc_src)
            pc_ref.generate('protein_coding')

            lncRNA_ref = EnsemblTranscriptsFile(self.lncRNA_src)
            lncRNA_ref.generate('lncRNA')

        elif self.src_type == 'gencode':
            pc_ref = self.pc_src
            lncRNA_ref = self.lncRNA_src

        repChrM_ref = RepChrM(self.repChrM_src)
        repChrM_ref.generate()

        with open(self.filename, 'wb') as out:
            for ref in [pc_ref, lncRNA_ref, repChrM_ref]:
                with open(ref.filename, 'rb') as f_in:
                    shutil.copyfileobj(f_in, out)

        return self.filename


@contextmanager
def cwd(path):
    origin_pwd = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(origin_pwd)


def generate(species, source, version, ref_dir):
    with cwd(ref_dir):
        species = species_list[species]

        if source == "ensembl":
            anno_file = rs.EnsemblAnnotation(species.name, version)
            genome_file = rs.EnsemblGenome(species.name, version)
            other_transcripts_files = (
                'ensembl',
                rs.EnsemblCDna(species.name, version),
                rs.EnsemblNCRna(species.name, version)
            )

        elif source == "gencode":

            if species.key == 'hsa':
                species_key = 'human'
            elif species.key == 'mmu':
                species_key = 'mouse'

            anno_file = rs.GencodeAnnotation(species_key, version)
            genome_file = rs.GencodeGenome(species_key, version)
            other_transcripts_files = (
                'gencode',
                rs.GencodeProteinCodingTranscripts(species_key, version),
                rs.GencodeLongNonCodingTranscripts(species_key, version)
            )

        elif source.startswith("ensembl_"):
            field = source.split('_')[1]
            anno_file = rs.EnsemblSisterAnnotation(field, species.name, version)
            genome_file = rs.EnsemblSisterGenome(field, species.name, version)
            other_transcripts_files = (
                'ensembl',
                rs.EnsemblSisterCDna(field, species.name, version),
                rs.EnsemblNCRna(field, species.name, version)
            )

        else:
            raise rs.SourceNotSupportError(source)

        mir_seq_file = rs.MiRBaseMiRNA(None, "21")
        mir_taret_file = rs.MiRTarBaseResource(None, "7.0")

        # download
        anno_file.download()
        genome_file.download()
        mir_seq_file.download()
        mir_taret_file.download()

        for file_ in other_transcripts_files[1:]:
            file_.download()

        # unzip
        genome_file.unzip()
        mir_seq_file.unzip()

        for file_ in other_transcripts_files[1:]:
            file_.unzip()

        # genref
        anno_ref = AnnoRef(anno_file.filename)
        anno_ref.generate()

        mir_ref = MirRef(mir_seq_file.filename)
        mir_ref.generate(species)

        mir_taret_ref = MiRTarBaseRef(mir_taret_file.filename)
        mir_taret_ref.generate(species)

        others_ref = OtherTranscripts(*other_transcripts_files, genome_file)
        others_ref.generate()

        # config
        info = {
            'species': species.key,
            'source': source,
            'version': anno_file.version
        }

        ref_files = {
            'anno_db': anno_ref.filename,
            'ref_file': genome_ref.filename,
            'mir_ref': mir_ref.filename,
            'mir_target': mir_taret_ref.filename,
            'other_transcripts': others_ref.filename
        }

        return info, ref_files
