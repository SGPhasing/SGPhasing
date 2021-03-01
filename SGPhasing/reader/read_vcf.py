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
  - get_pos_ref_alt
"""


def get_pos_ref_alt(input_vcf: str, ploidy: int) -> dict:
    """Get position ref base and alt base.

    Args:
        input_vcf (str): input vcf file path string.
        ploidy (int): ploidy number.

    Returns:
        pos_ref_alt (dict): position as key and [REF, ALT] as value.
    """
    min_freq = 0.5 / ploidy
    max_freq = 1 - min_freq
    pos_ref_alt = {}
    with open(input_vcf, 'r') as opened_vcf:
        for eachline in opened_vcf:
            if eachline[0] != '#':
                sp = eachline.strip().split()
                if len(sp[4]) == 1:
                    merge_sp = sp[9].split(':')
                    dp = int(merge_sp[2])
                    if all([min_freq < int(ph_dp) / dp < max_freq
                            for ph_dp in merge_sp[1].split(',')]):
                        if len(sp[3]) == 1:
                            pos_ref_alt.update({int(sp[1])-1: sp[3:5]})
                        else:
                            pos_ref_alt.update({int(sp[1]): [sp[3][1:], '-']})
    return pos_ref_alt
