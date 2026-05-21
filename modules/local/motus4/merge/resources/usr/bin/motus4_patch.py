#!/usr/bin/env python3
# mOTUs 4 does not support anymore the -db flag to pass a custom database
# the database path is hardcoded in the module, and this patch overrides it to
# point to the local work directory

import sys,pathlib
from unittest.mock import patch

workdir = pathlib.Path(__file__).resolve().parent

with (
    patch('motus.mutils.DEFAULT_MOTUS_MGDB_PARENT_LOCATION', new=workdir),
    patch('motus.mutils.DEFAULT_MOTUS_MGDB_LOCATION', new=workdir.joinpath('db_mOTU')),
    patch('motus.mutils.DEFAULT_MOTUS_MGDB_LOCATION_MARKER', new=workdir.joinpath('db_mOTU/db_mOTU.downloaded'))
):
    if __name__ == '__main__':
        sys.argv[0] = sys.argv[0].removesuffix('.exe')
        sys.exit(main())