# -*- coding: utf-8 -*-

"""Run this script with :code:`python3 -m bio2bel_drugbank`"""

import click

from bio2bel import build_cli
from .manager import Manager

main = build_cli(Manager)


@main.group()
def manage():
    pass


@manage.group()
def patents():
    pass


@patents.command()
@click.pass_obj
def ls(manager):
    """Lists patents"""
    click.echo_via_pager('\n'.join(
        f'{patent.patent_id}\t{patent.country}\t{"|".join(drug.name for drug in patent.drugs)}'
        for patent in manager.list_patents()
    ))


if __name__ == '__main__':
    main()
