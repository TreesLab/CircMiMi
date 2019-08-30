#! /usr/bin/env python
import click
from circmimi import Circmimi


@click.group()
def cli():
    pass


@cli.command()
@click.argument('circ_file')
@click.option('-a', '--annotation', 'anno_db')
@click.option('-g', '--genome', 'ref_file')
@click.option('-mir', '--mirna', 'mir_ref')
@click.option('-p', '--num_proc', default=1, type=click.INT)
@click.option('--header', 'header', flag_value=True, type=click.BOOL, default=True, hidden=True)
@click.option('--no-header', 'header', flag_value=False, type=click.BOOL)
def run(circ_file, anno_db, ref_file, mir_ref, num_proc, header):
    circmimi_result = Circmimi()
    circmimi_result.run(
        circ_file,
        anno_db,
        ref_file,
        mir_ref,
        num_proc=num_proc
    )

    result_table = circmimi_result.get_result_table()
    print(result_table.to_csv(sep='\t', index=False, header=header), end='')


if __name__ == "__main__":
    cli()
