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

    def create(self, ref_dir='.'):
        dat_file_from_ = rs.MiRBaseDat(self.species, self.from_)
        dat_file_to_ = rs.MiRBaseDat(self.species, self.to_)
        diff_file = rs.MiRBaseDiff(self.species, self.to_)

        dat_file_from_.download(ref_dir)
        dat_file_from_.rename_with_version()
        dat_file_to_.download(ref_dir)
        dat_file_to_.rename_with_version()
        diff_file.download(ref_dir)

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
            r'^([^\t]*)\t({}-[^\t]*)\t(.*)'.format(species),
            mature_diff_text,
            flags=re.M
        )

        return diff_of_species

    @staticmethod
    def _get_accessions_of_type(diff_type, diff_data):
        diff_data_of_type = filter(lambda data: diff_type in data[2], diff_data)
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

    def save(self, mapping_file):
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

    def update_row(self, row_data, col_key=0, inplace=False, remove_deleted=False):
        updated_row = list(row_data)
        updated_value = self.update(updated_row[col_key])

        if remove_deleted:
            if updated_value == '':
                return

        if inplace:
            updated_row[col_key] = updated_value
        else:
            updated_row.append(updated_value)

        return updated_row

    def update_file(self, in_file, out_file, col_key=0, inplace=False, remove_deleted=False):
        with open(in_file) as f_in, open(out_file, 'w') as f_out:
            for line in f_in:
                data = line.rstrip('\n').split('\t')
                updated_data = self.update_row(data, col_key, inplace, remove_deleted)
                if updated_data is not None:
                    print(*updated_data, sep='\t', file=f_out)
