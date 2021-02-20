# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.reader read sam, bam or cram file."""

from sys import exit

import pysam

from SGPhasing.sys_output import Output


def open_xam(input_xam: str) -> tuple:
    """Check input and open.

    Args:
        input_xam (str): input sam or bam file path str.

    Returns:
        xamfile (pysam.AlignmentFile): pysam open format.
        input_format (str): return input_xam format, sam or bam.
    """
    if input_xam.endswith('cram'):
        xamfile = pysam.AlignmentFile(input_xam, 'rc', check_sq=False)
        input_format = 'cram'
    elif input_xam.endswith('bam'):
        xamfile = pysam.AlignmentFile(input_xam, 'rb', check_sq=False)
        input_format = 'bam'
    elif input_xam.endswith('sam'):
        xamfile = pysam.AlignmentFile(input_xam, 'r', check_sq=False)
        input_format = 'sam'
    else:
        output = Output()
        output.error('input error: input format must be'
                     ' cram, bam or sam file.')
        exit()
    return xamfile, input_format


def check_index(input_xam: str, threads: int = 1) -> None:
    """Build index for bam if not exists.

    Args:
        input_xam (str): input sam or bam file path str.
        threads (int): threads using for pysam sort and index.
    """
    xamfile, input_format = open_xam(input_xam)
    try:
        xamfile.check_index()
        return input_xam
    except ValueError:
        output = Output()
        output.info(f'preparing samtools index for input {input_xam}.')
        pysam.index(input_xam, '-b', '-m', '17', '-@', str(threads))
        xamfile.close()
        return input_xam
    except AttributeError:
        output = Output()
        output.info(f'preparing samtools index for input {input_xam}.')
        if input_format == 'sam':
            input_bam = input_xam[:-3] + 'bam'
            pysam.sort('-o', input_bam, '--output-fmt', 'BAM',
                       '--threads', str(threads), input_xam)
            pysam.index(input_bam, '-b', '-m', '17', '-@', str(threads))
            xamfile.close()
            return input_bam
        else:
            output.error('input error: input format must be'
                         ' cram, bam or sam file.')
            exit()


def check_flag(read: pysam.AlignedSegment, flag: int) -> bool:
    """Check read if contain the flag.

    Args:
        read (pysam.AlignedSegment): input a read pysam.AlignedSegment.
        flag (int): input a flag.

    Returns:
        has_flag (bool): read has the flag or not.
    """
    if read.flag < flag:
        return False
    else:
        flag_bin, read_bin = bin(flag), bin(read.flag)
        for bin_id in range(1, len(flag_bin)-1):
            if flag_bin[-bin_id] == 1:
                if read_bin[-bin_id] != 1:
                    return False
        else:
            return True


def read_to_fastq(read: pysam.AlignedSegment) -> str:
    """Turn a pysam read to fastq 4 lines string.

    Args:
        read (pysam.AlignedSegment): input a read pysam.AlignedSegment.

    Returns:
        read (str): read has the flag or not.
    """
    return '\n'.join(['@'+read.query_name, read.query_sequence,
                      '+', pysam.array_to_qualitystring(read.query_qualities)])
