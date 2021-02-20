# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""Represent a full help argument parser and execute.

What's here:

Loads the relevant script modules and executes the script.
----------------------------------------------------------

Classes:
  - Index
"""

from gc import collect
from logging import getLogger
from pathlib import Path
# from subprocess import Popen

from SGPhasing.processor.collapse import collapse_isoforms_by_sam
from SGPhasing.reader import read_bed, read_fastx, read_xam
from SGPhasing.writer import write_fastx, write_xam
from SGPhasing.sys_output import Output

logger = getLogger(__name__)  # pylint: disable=invalid-name


class Index(object):
    """The Index process.

    Attributes:
      - args: Arguments.
      - output: Output info, warning and error.
    """

    def __init__(self, arguments) -> None:
        """Initialize Index.

        Args:
          - arguments: Arguments.
        """
        self.args = arguments
        self.output = Output()
        if not self.args.tmp:
            self.args.tmp = self.args.input + '.sgphasing.tmp'
        self.output.info(
            f'Initializing {self.__class__.__name__}: (args: {arguments}')
        logger.debug(
            f'Initializing {self.__class__.__name__}: (args: {arguments}')

    def prepare(self) -> None:
        """Check genome index."""
        self.tmp_floder_path = Path(self.args.tmp)
        if not self.tmp_floder_path.is_dir():
            self.tmp_floder_path.mkdir()
            self.output.info(f'creating temporary folder at {self.args.tmp}.')
        read_fastx.check_index(self.args.ref, self.args.threads)
        self.args.input = read_xam.check_index(self.args.input)

    def check_bed(self) -> None:
        """Check input limitation bed file."""
        if self.args.bed:
            self.limit_region_dict = read_bed.open_bed(self.args.bed)
        else:
            self.limit_region_dict = {}

    def get_multimapped_reads(self) -> None:
        """Get multiply mapped reads from input bam."""
        self.multimapped_reads_set = set()
        self.opened_xam, input_format = read_xam.open_xam(self.args.input)
        for read in self.opened_xam.fetch():
            if read.is_secondary:
                self.multimapped_reads_set.add(read.query_name)
        if self.limit_region_dict:
            region_reads_set = set()
            for chrom, region_list in self.limit_region_dict.items():
                for start, end in region_list:
                    for read in self.opened_xam.fetch(chrom, start, end):
                        region_reads_set.add(read.query_name)
            self.multimapped_reads_set &= region_reads_set
            del region_reads_set
            collect()

    def get_primary_region(self) -> None:
        """Get primary region from multiply mapped reads."""
        self.primary_sam_path = self.tmp_floder_path / 'primary_reference.sam'
        (self.primary_region,
            self.primary_reads_set) = write_xam.write_partial_sam(
                self.opened_xam, str(self.primary_sam_path),
                self.limit_region_dict, self.multimapped_reads_set)
        del self.multimapped_reads_set
        collect()

        self.opened_fastx, self.fastx_format = read_fastx.open_fastx(
            self.args.fastx)
        self.primary_fastx_path = (self.tmp_floder_path /
                                   ('primary.' + self.fastx_format))
        write_fastx.write_partial_fastx(
            self.opened_fastx, self.primary_fastx_path.open('w'),
            self.fastx_format, self.primary_reads_set)
        self.opened_fastx.close()

    def collapse_primary_sam(self):
        self.primary_gff, self.primary_group = collapse_isoforms_by_sam(
            str(self.primary_sam_path),
            str(self.primary_fastx_path),
            self.fastx_format == 'fastq',
            self.tmp_floder_path)

    def process(self) -> None:
        """Call the index object."""
        logger.debug('Starting index Process')
        self.prepare()
        self.check_bed()
        self.get_multimapped_reads()
        self.get_primary_region()
        self.collapse_primary_sam()
        logger.debug('Completed index Process')
