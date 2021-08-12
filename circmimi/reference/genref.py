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
from circmimi.reference.mirbase import MatureMiRNAUpdater


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
        df = pd.read_excel(self.src_name, engine='openpyxl')

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


class Files:
    def __init__(self, files):
        self.files = files

    def download(self):
        for file_ in self.files:
            file_.download()

    def unzip(self):
        for file_ in self.files:
            file_.unzip()

    def __getitem__(self, n):
        return self.files[n]


class EnsemblFiles(Files):
    source = "ensembl"


class GencodeFiles(Files):
    source = "gencode"


class MirTargetRef:
    def __init__(self, ref_files, ref_names, species):
        self.ref_files = ref_files
        self.ref_names = ref_names
        self.species = species
        self.filename = "mir_target.{}.tsv".format(self.species.key)

        self.remove_unavailable_files()

    def generate(self):
        merged_df = pd.DataFrame([], columns=['mirna', 'target_gene'])
        for ref_file, ref_name in zip(self.ref_files, self.ref_names):
            ref_df = pd.read_csv(ref_file, sep='\t', dtype='object')
            ref_df = ref_df.pipe(self.add_ref_name, ref_name).pipe(self.add_ref_col, ref_name)
            merged_df = merged_df.merge(ref_df, on=['mirna', 'target_gene'], how="outer")

        for ref_name in self.ref_names:
            merged_df[ref_name].fillna('0', inplace=True)

        # promote the refname columns
        col_names = list(merged_df.columns)
        col_names = col_names[:2] + self.promote_items(col_names[2:], self.ref_names)
        merged_df = merged_df[col_names]

        merged_df.to_csv(self.filename, sep='\t', index=False)

    @staticmethod
    def add_ref_name(ref_df, ref_name, sep='__'):
        data_cols = ref_df.columns[2:]
        new_col_names = dict(map(lambda col: (col, "{}{}{}".format(ref_name, sep, col)), data_cols))
        return ref_df.rename(new_col_names, axis=1)

    @staticmethod
    def add_ref_col(ref_df, ref_name):
        return ref_df.assign(**{ref_name: '1'})

    @staticmethod
    def promote_items(all_items, to_be_promoted):
        result_items = all_items.copy()
        for item in to_be_promoted[::-1]:
            if item in result_items:
                item_idx = result_items.index(item)
                result_items = [result_items[item_idx]] + result_items[:item_idx] + result_items[(item_idx + 1):]
        return result_items

    def remove_unavailable_files(self):
        ref_files = []
        ref_names = []
        for ref_file, ref_name in zip(self.ref_files, self.ref_names):
            if ref_file:
                ref_files.append(ref_file)
                ref_names.append(ref_name)

        self.ref_files = ref_files
        self.ref_names = ref_names


def generate(species, source, version, ref_dir):
    with cwd(ref_dir):
        species = species_list[species]

        if source == "ensembl":
            anno_file = rs.EnsemblAnnotation(species.name, version)
            genome_file = rs.EnsemblGenome(species.name, version)
            other_transcripts_files = EnsemblFiles(
                [
                    rs.EnsemblCDna(species.name, version),
                    rs.EnsemblNCRna(species.name, version)
                ]
            )

        elif source == "gencode":

            if species.key == 'hsa':
                species_key = 'human'
            elif species.key == 'mmu':
                species_key = 'mouse'

            anno_file = rs.GencodeAnnotation(species_key, version)
            genome_file = rs.GencodeGenome(species_key, version)
            other_transcripts_files = GencodeFiles(
                [
                    rs.GencodeProteinCodingTranscripts(species_key, version),
                    rs.GencodeLongNonCodingTranscripts(species_key, version)
                ]
            )

        elif source.startswith("ensembl_"):
            field = source.split('_')[1]
            anno_file = rs.EnsemblSisterAnnotation(field, species.name, version)
            genome_file = rs.EnsemblSisterGenome(field, species.name, version)
            other_transcripts_files = EnsemblFiles(
                [
                    rs.EnsemblSisterCDna(field, species.name, version),
                    rs.EnsemblSisterNCRna(field, species.name, version)
                ]
            )

        else:
            raise rs.SourceNotSupportError(source)

        mir_seq_file = rs.MiRBaseMiRNA(None, "22")
        mir_target_files = Files(
            [
                rs.MiRTarBaseResource(None, "7.0"),
                rs.MiRDBData(species.key, "6.0"),
                rs.EncoriMiRNATargetData(species.key)
            ]
        )
        ENCORI_RBP_files = Files(
            [
                rs.EncoriRBPData(species.key, source, version, only_AGO=True),
                # rs.EncoriRBPData(species.key, source, version),
                # rs.EncoriRBPTargetData(species.key)
            ]
        )

        # download
        anno_file.download()
        genome_file.download()
        mir_seq_file.download()

        mir_target_files.download()
        other_transcripts_files.download()
        ENCORI_RBP_files.download()

        # unzip
        genome_file.unzip()
        mir_seq_file.unzip()

        mir_target_files.unzip()
        other_transcripts_files.unzip()
        ENCORI_RBP_files.unzip()

        # genref
        anno_ref = AnnoRef(anno_file.filename)
        anno_ref.generate()

        mir_ref = MirRef(mir_seq_file.filename)
        mir_ref.generate(species.key)

        # miRTarBase
        miRTarBase_ref = MiRTarBaseRef(mir_target_files[0].filename)
        miRTarBase_ref.generate(species)

        # miRNAs updater
        updater = MatureMiRNAUpdater("21", "22", species.key)
        updater.create()

        updated_miRTarBase_ref_filename = "miRTarBase.{}.miRBase_v22.tsv".format(species.key)
        updater.update_file(
            miRTarBase_ref.filename,
            updated_miRTarBase_ref_filename,
            col_key=0,
            inplace=True,
            remove_deleted=True
        )

        mir_target_refs = [
            updated_miRTarBase_ref_filename,
            mir_target_files[1].filename,
            mir_target_files[2].filename
        ]
        mir_target_ref = MirTargetRef(
            mir_target_refs,
            [
                "miRTarBase",
                "miRDB",
                "ENCORI"
            ],
            species
        )
        mir_target_ref.generate()

        others_ref = OtherTranscripts(
            other_transcripts_files.source,
            *other_transcripts_files.files,
            genome_file
        )
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
            'other_transcripts': others_ref.filename,
            'AGO_data': ENCORI_RBP_files[0].filename,
            # 'RBP_data': ENCORI_RBP_files[1].filename,
            # 'RBP_target': ENCORI_RBP_files[2].filename
        }

        return info, ref_files
