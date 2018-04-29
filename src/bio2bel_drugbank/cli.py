# -*- coding: utf-8 -*-

"""Run this script with :code:`python3 -m bio2bel_drugbank`."""

import sys

import click

from .manager import Manager

main = Manager.get_cli()


@main.group()
def manage():
    """Manage the database."""


@manage.group()
def patents():
    """Manage patents."""


@patents.command()
@click.pass_obj
def ls(manager):
    """List patents."""
    click.echo_via_pager('\n'.join(
        f'{patent.patent_id}\t{patent.country}\t{"|".join(drug.name for drug in patent.drugs)}'
        for patent in manager.list_patents()
    ))


@manage.group()
def xrefs():
    """Manage cross-references."""


@xrefs.command()
@click.pass_obj
def summarize(manager):
    """List cross-reference types."""
    for resource, count in sorted(manager.summarize_xrefs()):
        click.echo(f'{resource}: {count}')


@xrefs.command()
@click.argument('resource')
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
@click.pass_obj
def build_mapping(manager, resource, output):
    """Builds a mapping file (TSV) for the given resource."""
    click.echo(f'drugbank\t{resource}', file=output)
    for xref in manager.get_xrefs_by_resource(resource):
        click.echo(f'{xref.drug.drugbank_id}\t{xref.identifier}', file=output)


if __name__ == '__main__':
    main()
