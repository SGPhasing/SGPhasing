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

from logging import getLogger
from pathlib import Path
# from subprocess import Popen

from SGPhasing.reader import read_bed, read_fasta, read_xam
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
        self.output.info(
            f'Initializing {self.__class__.__name__}: (args: {arguments}')
        logger.debug(
            f'Initializing {self.__class__.__name__}: (args: {arguments}')

    def prepare(self) -> None:
        """Check genome index."""
        if not self.args.tmp:
            self.args.tmp = self.args.input + '.sgphasing.tmp'
        self.tmp_floder = Path(self.args.tmp)
        if not self.tmp_floder.is_dir():
            self.tmp_floder.mkdir()
            self.output.info(f'creating temporary folder at {self.args.tmp}.')
        # minimap2_index = Path(self.args.ref + '.mmi')
        # if not minimap2_index.exists():
        #     with Popen(['minimap2', '-k', '17', '-I', '18G', '-x', 'splice',
        #                 '-d', f'{self.args.ref}.mmi', self.args.ref]):
        read_fasta.check_index(self.args.ref, self.args.threads)
        self.output.info('preparing genome index for minimap2.')
        read_xam.check_index(self.args.input)

    def check_bed(self) -> None:
        """Check input limitation bed file."""
        if self.args.bed:
            self.limit_region = read_bed.open_bed(self.args.bed)
        else:
            self.limit_region = {}

    def get_multimapped_reads(self) -> None:
        """Get multiply mapped reads from input bam."""
        self.multimapped_reads_set = set()
        self.opened_xam, input_format = read_xam.open_xam(self.args.input)
        for read in self.opened_xam.fetch():
            if read.is_secondary:
                self.multimapped_reads_set.add(read.query_name)
        if self.limit_region:
            region_reads_set = set()
            for chrom, region_list in self.limit_region.items():
                for start, end in region_list:
                    for read in self.opened_xam.fetch(chrom, start, end):
                        region_reads_set.add(read.query_name)
            self.multimapped_reads_set &= region_reads_set
            del region_reads_set

    def get_primary_region(self) -> None:
        """Get primary region from multiply mapped reads."""
        if self.limit_region:
            for chrom, region_list in self.limit_region.items():
                for start, end in region_list:
                    for read in self.opened_xam.fetch(chrom, start, end):
                        if read.query_name in self.multimapped_reads_set:
                            if not read.is_secondary:
                                read.reference_start, read.reference_end


    def process(self) -> None:
        """Call the index object."""
        logger.debug('Starting index Process')
        self.prepare()
        self.check_bed()
        self.get_multimapped_reads()
        logger.debug('Completed index Process')
