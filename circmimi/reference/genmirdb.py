import pandas as pd
from functools import partial
from circmimi.reference.utils import cwd, open_file
from circmimi.reference.resource import MiRDBResource, Gene2AccessionResource
from circmimi.reference.species import species_list


def get_needed_data(filename, data_filter):
    with open_file(filename) as f_in:
        for line in data_filter(f_in):
            yield line.rstrip('\n').split('\t')


def mirdb_filter(data, species_key):
    for line in data:
        if line.startswith(species_key):
            yield line


def gene_accession_filter(data, species_tax_id):
    species_tax_id_txt = "{}\t".format(species_tax_id)

    title = data.readline()
    yield title

    for line in data:
        if line.startswith(species_tax_id_txt):
            yield line


def generate(species, version, ref_dir, out_file, show_accession=False):
    species = species_list[species]

    with cwd(ref_dir):
        mirdb_file = MiRDBResource(None, version)
        gene2accession = Gene2AccessionResource()

        # download
        mirdb_file.download()
        gene2accession.download()

        # basic filter for the "species"
        ## mirdb
        mirdb_data = get_needed_data(
            mirdb_file.filename,
            data_filter=partial(mirdb_filter, species_key=species.key)
        )

        mirdb_df = pd.DataFrame(
            mirdb_data,
            columns=['mirna', 'accession', 'targeting_score']
        )

        ## gene2accesion
        gene_accession_data = get_needed_data(
            gene2accession.filename,
            data_filter=partial(
                gene_accession_filter,
                species_tax_id=species.tax_id
            )
        )

        gene_accession_data_title = next(gene_accession_data)

        gene_accession_df = pd.DataFrame(
            gene_accession_data,
            columns=gene_accession_data_title
        ).pipe(
            lambda df: df[['RNA_nucleotide_accession.version', 'Symbol', 'GeneID']]
        ).drop_duplicates(
        ).pipe(
            lambda df: df[df['RNA_nucleotide_accession.version'] != '-']
        )

        accession_no_version_df = gene_accession_df.apply(
            lambda s: pd.Series(s['RNA_nucleotide_accession.version'].split('.')[0]),
            axis=1
        )

        gene_accession_df_2 = gene_accession_df.assign(
            accession=accession_no_version_df
        ).pipe(
            lambda df: df[['accession', 'Symbol']]
        ).rename(
            {
                'Symbol': 'target_gene'
            },
            axis=1
        )

        # merge
        mirdb_df_with_symbol = mirdb_df.merge(
            gene_accession_df_2,
            on='accession',
            how='inner'
        ).pipe(
            lambda df: df[['mirna', 'target_gene', 'targeting_score', 'accession']]
        )

        if not show_accession:
            mirdb_df_with_symbol = mirdb_df_with_symbol.drop(
                'accession',
                axis=1
            ).groupby(
                [
                    'mirna',
                    'target_gene'
                ]
            ).agg('max').reset_index()

        mirdb_df_with_symbol.to_csv(out_file, sep='\t', index=None)
