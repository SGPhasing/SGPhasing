# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.reader read bam, cram or sam file.

Functions:
  - open_xam
  - check_index
  - extract_read_matrix
  - remove_blank
  - check_flag
  - read_to_fastq
"""

from sys import exit

import pysam

from SGPhasing.sys_output import Output


def open_xam(input_xam: str) -> tuple:
    """Check input and open.

    Args:
        input_xam (str): input sam or bam file path string.

    Returns:
        xamfile (AlignmentFile): pysam opened bam/cram/sam file handle.
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
        output.error('Input error: input format must be'
                     ' bam, cram or sam file.')
        exit()
    return xamfile, input_format


def check_index(input_xam: str, threads: int = 1) -> str:
    """Build index for bam if not exists.

    Args:
        input_xam (str): input sam or bam file path string.
        threads (int): threads using for pysam sort and index, default 1.
    """
    xamfile, input_format = open_xam(input_xam)
    try:
        xamfile.check_index()
        return input_xam
    except ValueError:
        output = Output()
        output.info(f'Preparing samtools index for input {input_xam}')
        pysam.index(input_xam, '-b', '-m', '17', '-@', str(threads))
        xamfile.close()
        return input_xam
    except AttributeError:
        output = Output()
        output.info(f'Preparing samtools index for input {input_xam}')
        if input_format == 'sam':
            input_bam = input_xam[:-3] + 'bam'
            pysam.sort('-o', input_bam, '--output-fmt', 'BAM',
                       '--threads', str(threads), input_xam)
            pysam.index(input_bam, '-b', '-m', '17', '-@', str(threads))
            xamfile.close()
            return input_bam
        else:
            output.error('Input error: input format must be'
                         ' bam, cram or sam file.')
            exit()


def extract_read_matrix(input_bam: str, positions_list: list) -> tuple:
    """Extract reads base matrix.

    Args:
        input_bam (str): input bam file path string.
        positions_list (list): position list for each base.

    Returns:
        reads_id_list (list): extracted read_id list for each read.
        read_base_matrix (list): extracted bases matrix for each read.
    """
    reads_bases_matrix, reads_id_list = [], []
    opened_bam = pysam.AlignmentFile(input_bam, 'rb')
    bases_num = len(positions_list)
    for read_index, read in enumerate(opened_bam.fetch()):
        read_id = read.query_name
        reads_id_list.append(read_id)
        reads_bases_matrix.append(['' for base_id in range(bases_num)])
    base_id = 0
    for pileupcolumn in opened_bam.pileup():
        if pileupcolumn.pos in positions_list:
            for pileupread in pileupcolumn.pileups:
                read_index = reads_id_list.index(
                    pileupread.alignment.query_name)
                if pileupread.is_del:
                    reads_bases_matrix[read_index][base_id] = '-'
                else:
                    reads_bases_matrix[read_index][base_id] = (
                        pileupread.alignment.query_sequence[
                            pileupread.query_position])
            base_id += 1
    opened_bam.close()
    return reads_id_list, reads_bases_matrix


def remove_blank(reads_id_list: list, reads_bases_matrix: list) -> tuple:
    """Remove blank reads from reads_bases_matrix.

    Args:
        reads_id_list (list): extracted read_id list for each read.
        read_base_matrix (list): extracted bases matrix for each read.

    Returns:
        new_reads_id_list (list): blank-removed read_id list for each read.
        new_reads_bases_matrix (list): blank-removed bases matrix
                                       for each read.
    """
    new_reads_id_list, new_reads_bases_matrix = [], []
    for read_id, read_bases_matrix in zip(reads_id_list, reads_bases_matrix):
        if not all([base == '' for base in read_bases_matrix]):
            new_reads_id_list.append(read_id)
            new_reads_bases_matrix.append(read_bases_matrix)
    return new_reads_id_list, new_reads_bases_matrix


def check_flag(read: pysam.AlignedSegment, flag: int) -> bool:
    """Check read if contain the flag.

    Args:
        read (AlignedSegment): read in pysam.AlignedSegment type.
        flag (int): flag.

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
        read (AlignedSegment): read in pysam.AlignedSegment type.

    Returns:
        read (str): read has the flag or not.
    """
    return '\n'.join(['@'+read.query_name, read.query_sequence,
                      '+', pysam.array_to_qualitystring(read.query_qualities)])
