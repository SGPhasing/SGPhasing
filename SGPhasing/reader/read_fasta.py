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

import mappy as mp


def check_index(input_ref: str, threads: int = 1) -> None:
    opened_aligner = mp.Aligner(input_ref, preset='splice', k=17, best_n=100,
                                n_threads=threads, fn_idx_out=input_ref+'.mmi')
