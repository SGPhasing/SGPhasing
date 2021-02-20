# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.writer write xam file."""

from pysam import AlignmentFile

from SGPhasing.reader.read_bed import merge_region


def write_partial_sam(opened_input_xam: AlignmentFile,
                      output_sam: str,
                      limit_region_dict: dict,
                      limit_reads_set: set) -> tuple:
    """Write sam for limit region and reads.

    Args:
        opened_input_xam (pysam.AlignmentFile): pysam open format.
        output_sam (str): output sam path str.
        limit_region_dict (dict): chrom as key and region list as value.
        limit_reads_set (set): limited reads id set.

    Returns:
        chr_region (dict): chrom as key and region list as value.
        output_reads_set (set): output reads id set.
    """
    opened_output_sam = AlignmentFile(
        output_sam, 'w', template=opened_input_xam)
    chr_region, output_reads_set = {}, set()
    if limit_region_dict:
        for chrom, region_list in limit_region_dict.items():
            for start, end in region_list:
                for read in opened_input_xam.fetch(chrom, start, end):
                    if (read.query_name in limit_reads_set and
                            not read.is_supplementary):
                        opened_output_sam.write(read)
                        output_reads_set.add(read.query_name)
                        chr_region.setdefault(chrom, []).append(
                            (read.reference_start, read.reference_end))
    else:
        for read in opened_input_xam.fetch():
            if (read.query_name in limit_reads_set and
                    not read.is_supplementary):
                opened_output_sam.write(read)
                output_reads_set.add(read.query_name)
                chr_region.setdefault(read.re, []).append(
                    (read.reference_start, read.reference_end))
    opened_output_sam.close()
    return merge_region(chr_region), output_reads_set
