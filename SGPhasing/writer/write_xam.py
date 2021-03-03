# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.writer write sam file.

Functions:
  - write_partial_sam
  - sam_to_bam
"""

import pysam

from SGPhasing.reader.read_bed import merge_region


def write_partial_sam(opened_input_xam: pysam.AlignmentFile,
                      output_sam: str,
                      limit_region_dict: dict,
                      limit_reads_set: set) -> tuple:
    """Write sam for limit region and reads.

    Args:
        opened_input_xam (AlignmentFile): input pysam opened
                                          bam/sam file handle.
        output_sam (str): output sam file path string.
        limit_region_dict (dict): chrom as key and region list as value.
        limit_reads_set (set): limited reads id set.

    Returns:
        chr_region (dict): chrom as key and region list as value.
        output_reads_set (set): output reads id set.
    """
    opened_output_sam = pysam.AlignmentFile(
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


def sam_to_bam(input_sam: str,
               reference: str,
               output_bam: str,
               threads: int = 1) -> None:
    """Sort sam and save to bam.

    Args:
        input_sam (str): input sam file path string.
        reference (str): input reference fasta file path string.
        output_bam (str): output bam file path string.
        threads (int): threads using for pysam.sort, default = 1.
    """
    pysam.sort('-o', output_bam, '--output-fmt', 'BAM',
               '--reference', reference, '--threads', str(threads), input_sam)
