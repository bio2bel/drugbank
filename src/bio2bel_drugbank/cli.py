# -*- coding: utf-8 -*-

"""Run this script with :code:`python3 -m bio2bel_drugbank`"""

from bio2bel import build_cli
from .manager import Manager

main = build_cli(Manager)

if __name__ == '__main__':
    main()
