#! /usr/bin/env python
import click
import os


@click.group()
def cli():
    pass


@cli.command(help="Main pipeline.")
@click.argument('circ_file')
@click.option('-r', '--ref', 'ref_dir', required=True, help="Assign the path of ref_dir.", type=click.Path())
@click.option('-o', '--out', 'output_dir', hidden=True) # Not implemented yet.
@click.option('-p', '--num_proc', default=1, type=click.INT, help="The number of processes.")
@click.option('--header', 'header', flag_value=True, type=click.BOOL,
              default=True, hidden=True)
@click.option('--no-header', 'header', flag_value=False, type=click.BOOL)
def run(circ_file, ref_dir, output_dir, num_proc, header):
    from circmimi.config import get_refs
    anno_db, ref_file, mir_ref, mir_target = get_refs(ref_dir)

    from circmimi.circmimi import Circmimi

    circmimi_result = Circmimi()
    circmimi_result.run(
        circ_file,
        anno_db,
        ref_file,
        mir_ref,
        mir_target,
        num_proc=num_proc
    )

    result_table = circmimi_result.get_result_table()
    print(result_table.to_csv(sep='\t', index=False, header=header), end='')


@cli.command()
@click.option('--species', 'species')
@click.option('--source', 'source')
@click.option('--gencode', 'source', flag_value='gencode', type=click.STRING)
@click.option('--ensembl', 'source', flag_value='ensembl', type=click.STRING)
@click.option('--version', 'version')
@click.option('--init', 'init', is_flag=True, help="Create an init template ref_dir.")
@click.argument('ref_dir')
def genref(species, source, version, ref_dir, init):
    os.makedirs(ref_dir, exist_ok=True)

    from circmimi.config import RefConfig
    config = RefConfig()

    if not init:
        # TO BE DONE!!
        print("Not implemented yet!")
        print("species: {}\nsource: {}\nversion: {}\nref_dir: {}".format(
            species, source, version, ref_dir))

    config.write(ref_dir)


@cli.command()
@click.argument('gtf_path')
@click.argument('db_path', metavar='OUT_PATH')
def gendb(gtf_path, db_path):
    from circmimi.reference import generate_db

    generate_db(gtf_path, db_path)


if __name__ == "__main__":
    cli()
