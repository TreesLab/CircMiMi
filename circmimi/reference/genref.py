import re
import os
import pandas as pd
import csv
import shutil
from circmimi.reference import gendb
from circmimi.reference.species import species_list
from circmimi.seq import parse_fasta
from circmimi.reference import resource as rs
from circmimi.reference.utils import cwd


class RefFile:
    def __init__(self, src_name):
        self.src_name = src_name
        self.filename = None

    def generate(self):
        return self.filename


class AnnoRef(RefFile):
    def generate(self):
        self.filename = re.sub(r'\.gtf(?:\.gz)?$', '.db', self.src_name)
        gendb.generate(self.src_name, self.filename)
        return self.filename


class MirRef(RefFile):
    def generate(self, species_key):
        self.filename = re.sub(r'\.fa$', '.{}.fa'.format(species_key), self.src_name)

        with open(self.src_name) as fa_in:
            fa_txt = fa_in.read()

        with open(self.filename, 'w') as out:
            for fa_id, fa_seq in parse_fasta(fa_txt):
                m = re.search(r'^{}-'.format(species_key), fa_id)
                if m:
                    print(">{}\n{}".format(fa_id, fa_seq), file=out)

        return self.filename


class MiRTarBaseRef(RefFile):
    def generate(self, species):
        df = pd.read_excel(self.src_name)

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
        self.filename = re.sub(r'\.fa$', '.{}.fa'.format(biotype), self.src_name)

        with open(self.src_name) as fa_in:
            fa_txt = fa_in.read()

        with open(self.filename, 'w') as out:
            for fa_id, fa_seq in parse_fasta(fa_txt):
                m = re.search(r'{}'.format(biotype), fa_id)
                if m:
                    print(">{}\n{}".format(fa_id, fa_seq), file=out)

        return self.filename


class RepChrM(RefFile):
    def generate(self):
        file_path = os.path.dirname(self.src_name)
        self.filename = os.path.join(file_path, "RepChrM.fa")

        with open(self.src_name) as fa_in:
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
            pc_ref = EnsemblTranscriptsFile(self.pc_src.filename)
            pc_ref.generate('protein_coding')

            lncRNA_ref = EnsemblTranscriptsFile(self.lncRNA_src.filename)
            lncRNA_ref.generate('lncRNA')

        elif self.src_type == 'gencode':
            pc_ref = self.pc_src
            lncRNA_ref = self.lncRNA_src

        repChrM_ref = RepChrM(self.repChrM_src.filename)
        repChrM_ref.generate()

        with open(self.filename, 'wb') as out:
            for ref in [pc_ref, lncRNA_ref, repChrM_ref]:
                with open(ref.filename, 'rb') as f_in:
                    shutil.copyfileobj(f_in, out)

        return self.filename


def generate(species, source, version, ref_dir, use_miRDB=False):
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
                rs.EnsemblSisterNCRna(field, species.name, version)
            )

        else:
            raise rs.SourceNotSupportError(source)

        if use_miRDB:
            mir_seq_file = rs.MiRBaseMiRNA(None, "22")
            mir_target_file = rs.MiRDBData(species.key, "6.0")
        else:
            mir_seq_file = rs.MiRBaseMiRNA(None, "21")
            mir_target_file = rs.MiRTarBaseResource(None, "7.0")

        # download
        anno_file.download()
        genome_file.download()
        mir_seq_file.download()
        mir_target_file.download()

        for file_ in other_transcripts_files[1:]:
            file_.download()

        # unzip
        genome_file.unzip()
        mir_seq_file.unzip()
        mir_target_file.unzip()

        for file_ in other_transcripts_files[1:]:
            file_.unzip()

        # genref
        anno_ref = AnnoRef(anno_file.filename)
        anno_ref.generate()

        mir_ref = MirRef(mir_seq_file.filename)
        mir_ref.generate(species.key)

        if use_miRDB:
            mir_target_ref = mir_target_file
        else:
            mir_target_ref = MiRTarBaseRef(mir_target_file.filename)
            mir_target_ref.generate(species)

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
            'ref_file': genome_file.filename,
            'mir_ref': mir_ref.filename,
            'mir_target': mir_target_ref.filename,
            'other_transcripts': others_ref.filename
        }

        return info, ref_files
