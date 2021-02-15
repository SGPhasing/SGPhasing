# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.reader read sam and bam file."""

from sys import exit

import pysam

from SGPhasing.sys_output import Output

output = Output()


def open_xam(input_xam: str) -> tuple:
    """Check input and open."""
    if input_xam.endswith('bam'):
        xamfile = pysam.AlignmentFile(input_xam, 'rb', check_sq=False)
        input_format = 'bam'
    elif input_xam.endswith('sam'):
        xamfile = pysam.AlignmentFile(input_xam, 'r', check_sq=False)
        input_format = 'sam'
    else:
        output.error('input error: input format must be bam or sam file.')
        exit()
    return xamfile, input_format


def check_index(input_xam: str, threads: int = 1) -> None:
    """Build index for bam if not exists."""
    xamfile, input_format = open_xam(input_xam)
    try:
        xamfile.check_index()
    except ValueError:
        input_bai = input_xam + '.bai'
        pysam.index(input_xam, input_bai)
    except AttributeError:
        if input_format == 'sam':
            input_bam = input_xam[:-3] + 'bam'
            pysam.sort('-o', input_bam,
                       '--threads', str(threads), input_xam)
            input_bai = input_bam + '.bai'
            pysam.index(input_bam, input_bai)
