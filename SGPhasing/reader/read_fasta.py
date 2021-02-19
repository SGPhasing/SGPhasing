# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.reader read fasta file."""

from pathlib import Path

import mappy as mp

from SGPhasing.sys_output import Output


def check_index(input_ref: str, threads: int = 1) -> None:
    """Build index for minimap2 if not exists.

    Args:
        input_ref (str): input reference fasta file path.
        threads (int): threads using for mappy to index.
    """
    index_file = input_ref + '.mmi'
    index_path = Path(index_file)
    if not index_path.exists():
        output = Output()
        output.info('preparing genome index for minimap2.')
        opened_aligner = mp.Aligner(input_ref, preset='splice',
                                    k=17, best_n=100,
                                    n_threads=threads, fn_idx_out=index_file)
