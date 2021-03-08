# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.reader read index file.

Functions:
  - read_index
  - array_to_seq
"""

from pathlib import Path

from h5py import File
import numpy as np

from SGPhasing.writer.write_fastx import write_fasta


def read_index(input_index: str, tmp_floder_path: Path) -> dict:
    """Read hdf5 index file.

    Args:
        input_index (str): input hdf5 index file path string.
        tmp_floder_path (Path): temporary folder Path.

    Returns:
        link_region_main_dict (dict): link id as key and
                                      main region dict as value.
    """
    link_region_main_dict = {}
    opened_index = File(input_index, 'r')
    for link_id in opened_index.keys():
        link_floder_path = tmp_floder_path / link_id
        link_fasta_path = (
            link_floder_path / 'expanded_primary.reference.fasta')
        if not link_floder_path.is_dir():
            link_floder_path.mkdir()
        for region_id in opened_index[link_id].keys():
            if region_id == 'reference':
                if not link_fasta_path.exists():
                    opened_link_fasta = link_fasta_path.open('w')
                    sequence = array_to_seq(
                        np.array(opened_index[link_id][region_id]))
                    write_fasta(opened_link_fasta, link_id, sequence)
                    opened_link_fasta.close()
            else:
                link_region_main_dict.setdefault(link_id, {}).update({
                    region_id: [
                        opened_index[link_id][region_id].attrs['chromosome'],
                        opened_index[link_id][region_id].attrs['start'],
                        opened_index[link_id][region_id].attrs['end'],
                        np.array(opened_index[link_id][region_id]['index'])]})
    opened_index.close()
    return link_region_main_dict


def array_to_seq(array: np.ndarray) -> str:
    """Convert numpy array to sequence.

    Args:
        array (ndarray): numpy array.

    Returns:
        (str): sequence string.
    """
    INT_BASE = ['A', 'C', 'G', 'T', 'M', 'R', 'W', 'S',
                'Y', 'K', 'V', 'H', 'D', 'B', 'N']
    return ''.join([INT_BASE[index] for index in array])
