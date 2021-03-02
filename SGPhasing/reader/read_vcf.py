# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.reader read vcf file.

Functions:
  - get_alt_positions
"""


def get_alt_positions(input_vcf: str, ploidy: int) -> list:
    """Get position ref base and alt base.

    Args:
        input_vcf (str): input vcf file path string.
        ploidy (int): ploidy number.

    Returns:
        positions_list (dict): positions list.
    """
    min_freq = 0.5 / ploidy
    max_freq = 1 - min_freq
    positions_list = []
    with open(input_vcf, 'r') as opened_vcf:
        for eachline in opened_vcf:
            if eachline[0] != '#':
                sp = eachline.strip().split()
                alt_sp = sp[4].split(',')
                if any([len(alt) == 1 for alt in alt_sp]):
                    merge_sp = sp[9].split(':')
                    dp = int(merge_sp[2])
                    if all([min_freq < int(ph_dp) / dp < max_freq
                            for ph_dp in merge_sp[1].split(',')[1:]]):
                        if len(sp[3]) == 1:
                            positions_list.append(int(sp[1])-1)
                        else:
                            positions_list.append(int(sp[1]))
    return positions_list
