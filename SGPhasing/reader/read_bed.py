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
    return merge_region(chr_region)


def merge_region(chr_region: dict) -> dict:
    """Merge overlap regions.

    Args:
        chr_region (dict): chrom as key and region list as value.

    Returns:
        chr_region (dict): chrom as key and region list as value.
    """
    return {chrom: fold_region(unfold_region(region_list))
            for chrom, region_list in chr_region.items()}


def unfold_region(region_list: list) -> set:
    """Turn region_list into position_set.

    Args:
        region_list (list): the elements are (start, end) tuple.

    Returns:
        position_set (set): set of positions.
    """
    position_set = set()
    for region_start, region_end in region_list:
        position_set |= set(range(region_start, region_end))
    return position_set


def fold_region(position_set: set) -> list:
    """Turn position_set into region_list.

    Args:
        position_set (set): set of positions.

    Returns:
        region_list (list): the elements are (start, end) tuple.
    """
    region_list = []
    sorted_pos = sorted(position_set)
    start = 0
    for pos_id in range(1, len(sorted_pos)):
        if sorted_pos[pos_id] - sorted_pos[pos_id-1] != 1:
            region_list.append((start, sorted_pos[pos_id-1]))
            start = sorted_pos[pos_id]
    return region_list
