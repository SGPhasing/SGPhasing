# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.writer write gff3 file.

Functions:
  - write_partial_gff
"""

from io import TextIOWrapper


def write_partial_gff(input_gff: str,
                      opened_output_gff: TextIOWrapper,
                      iso_id_list: list) -> None:
    """Write gff for limit isoforms.

    Args:
        input_gff (str): input gff file path string.
        opened_output_gff (TextIOWrapper): opened output gff handle.
        iso_id_list (list): limited isoforms id in list.
    """
    with open(input_gff, 'r') as opened_input_gff:
        for eachline in opened_input_gff:
            iso_id = eachline.strip().split()[11][1:-2]
            if iso_id in iso_id_list:
                opened_output_gff.write(eachline)
    opened_output_gff.close()
