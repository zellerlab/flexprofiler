#!/usr/bin/env python3
# mOTUs 4 does not support anymore the -db flag to pass a custom database
# the database path is hardcoded in the module, and this patch overrides it to
# point to the local work directory

# mOTUs 4 assigns an insert that aligns equally well to several genes of the
# same marker gene cluster to one of them at random, without a seed, so the
# profiles change slightly between runs. The seed makes them reproducible.

import sys,pathlib,random
from unittest.mock import patch
from motus.motus import main

random.seed(0)

workdir = pathlib.Path.cwd()

with (
    patch('motus.mutils.DEFAULT_MOTUS_MGDB_PARENT_LOCATION', new=workdir),
    patch('motus.mutils.DEFAULT_MOTUS_MGDB_LOCATION', new=workdir.joinpath('db_mOTU')),
    patch('motus.mutils.DEFAULT_MOTUS_MGDB_LOCATION_MARKER', new=workdir.joinpath('db_mOTU/db_mOTU.downloaded'))
):
    if __name__ == '__main__':
        sys.argv[0] = sys.argv[0].removesuffix('.exe')
        sys.exit(main())