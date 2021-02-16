# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.reader read bed file."""


def open_bed(input_bed: str) -> dict:
    """Read bed file.

    Args:
        input_bed (str): input bed file path.

    Returns:
        chr_region (dict): chrom as key and region list as value.
    """
    chr_region = {}
    with open(input_bed, 'r') as open_bed:
        for eachline in open_bed:
            if eachline[0] != '#':
                sp = eachline.strip().split()
                # chrom chromStart chromEnd name score strand
                chr_region.setdefault(sp[0], []).append(
                    (int(sp[1]), int(sp[2])))
    return chr_region
