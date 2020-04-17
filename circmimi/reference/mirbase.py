import re
from operator import itemgetter
from circmimi.reference import resource as rs
from circmimi.reference.utils import open_file


class MatureMiRNAUpdater:
    def __init__(self, from_, to_, species):
        self.from_ = from_
        self.to_ = to_
        self.species = species
        self.mapping_table = {}

    def create(self):
        dat_file_from_ = rs.MiRBaseDat(self.species, self.from_)
        dat_file_to_ = rs.MiRBaseDat(self.species, self.to_)
        diff_file = rs.MiRBaseDiff(self.species, self.to_)

        dat_file_from_.download()
        dat_file_from_.rename_with_version()
        dat_file_to_.download()
        dat_file_to_.rename_with_version()
        diff_file.download()

        dat_text_from_ = self._load_data(dat_file_from_.filename)
        dat_text_to_ = self._load_data(dat_file_to_.filename)
        diff_text = self._load_data(diff_file.filename)

        mature_diff_data = self._get_mature_diff(diff_text, self.species)

        renamed_accessions = self._get_accessions_of_type('NAME', mature_diff_data)
        renamed_miRNAs_from_ = self._get_IDs_by_accessions(dat_text_from_, renamed_accessions)
        renamed_miRNAs_to_ = self._get_IDs_by_accessions(dat_text_to_, renamed_accessions)
        self.mapping_table.update(dict(zip(renamed_miRNAs_from_, renamed_miRNAs_to_)))

        deleted_accessions = self._get_accessions_of_type('DELETE', mature_diff_data)
        deleted_miRNAs = list(self._get_IDs_by_accessions(dat_text_from_, deleted_accessions))
        self.mapping_table.update({mirna_id: '' for mirna_id in deleted_miRNAs})

    @staticmethod
    def _get_mature_diff(diff_text, species):
        mature_diff_text = re.search(
            r'#\n# Mature sequences start here\n#\n(.*)',
            diff_text,
            flags=re.DOTALL
        ).group(1)

        diff_of_species = re.findall(
            r'(.*)\t({}-.*)\t(.*)'.format(species),
            mature_diff_text
        )

        return diff_of_species

    @staticmethod
    def _get_accessions_of_type(diff_type, diff_data):
        diff_data_of_type = filter(lambda data: data[2] == diff_type, diff_data)
        accessions = list(map(itemgetter(0), diff_data_of_type))
        return accessions

    @staticmethod
    def _get_ID_by_accession(mirna_dat, accession):
        m = re.search(r'accession="{}"\n.*product="([^"]*)"'.format(accession), mirna_dat)
        if m:
            return m.group(1)
        else:
            return ''

    @classmethod
    def _get_IDs_by_accessions(cls, mirna_dat, accessions):
        for accession in accessions:
            yield cls._get_ID_by_accession(mirna_dat, accession)

    @staticmethod
    def _load_data(filename):
        with open_file(filename) as f_in:
            data_text = f_in.read()

        return data_text

    def save(self):
        mapping_file = 'miRNA.maps.{}_to_{}.tsv'.format(self.from_, self.to_)
        with open(mapping_file, 'w') as out:
            for k, v in self.mapping_table.items():
                print(k, v, sep='\t', file=out)

    def load_maps(self, mapping_file):
        with open(mapping_file) as f_in:
            for line in f_in:
                data = line.rstrip('\n').split('\t')
                self.mapping_table[data[0]] = data[1]

    def update(self, mirna):
        return self.mapping_table.get(mirna, mirna)
