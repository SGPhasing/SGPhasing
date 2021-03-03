# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""Index with alleles of high quality full-length transcriptome sequences.

What's here:

Index with alleles the process.
-------------------------------

Classes:
  - Index
"""

from gc import collect
from logging import getLogger
from multiprocessing import Pool
from pathlib import Path
import sys

from SGPhasing.processor.collapse import collapse_isoforms_by_sam
from SGPhasing.processor.gff_to_fasta import gff_to_fasta
from SGPhasing.processor.minimap2 import splice_mapper
from SGPhasing.processor.process_each_link import process_each_link
from SGPhasing.processor.sam_to_gff import sam_to_gff
from SGPhasing.reader import read_bed, read_fastx, read_gff, read_xam
from SGPhasing.writer import write_fastx, write_gff, write_xam
from SGPhasing.sys_output import Output

logger = getLogger(__name__)  # pylint: disable=invalid-name


class Index(object):
    """The Index process.

    Attributes:
        args: Arguments.
        output: Output info, warning and error.
    """

    def __init__(self, arguments) -> None:
        """Initialize Index.

        Args:
            arguments: arguments.
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
            self.output.info(f'Creating temporary folder at {self.args.tmp}')
        self.opened_log_file = (
            self.tmp_floder_path / 'sgphasing.log').open('w')
        read_fastx.check_index(self.args.reference, self.args.threads)
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
        self.primary_region, self.primary_reads_set = (
            write_xam.write_partial_sam(
                self.opened_xam, str(self.primary_sam_path),
                self.limit_region_dict, self.multimapped_reads_set))
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

    def collapse_primary_sam(self) -> None:
        """Collapse primary isoforms by sam."""
        self.output.info('Running cDNA_cupcake collapse_isoforms_by_sam')
        if not self.args.verbose:
            sys.stdout = self.opened_log_file
            sys.stderr = self.opened_log_file
            self.opened_log_file.write(
                'cDNA_cupcake collapse_isoforms_by_sam info:\n')
        self.primary_gff, self.most_iso_id_list = collapse_isoforms_by_sam(
            str(self.primary_sam_path), str(self.primary_fastx_path),
            self.fastx_format == 'fastq', self.tmp_floder_path)
        if not self.args.verbose:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__

    def get_most_supported_fasta(self) -> None:
        """Get most supported isoforms gff."""
        self.opened_most_gff = (self.tmp_floder_path /
                                'primary_reference.most.gff').open('w')
        write_gff.write_partial_gff(
            self.primary_gff, self.opened_most_gff, self.most_iso_id_list)
        self.primary_fasta_path = (self.tmp_floder_path /
                                   'primary_reference.most.fasta')
        self.output.info('Running cufflinks gffread')
        if self.args.verbose:
            self.output.info(gff_to_fasta(
                self.opened_most_gff.name, self.args.reference,
                str(self.primary_fasta_path)))
        else:
            self.opened_log_file.write('cufflinks gffread info:\n')
            self.opened_log_file.write(gff_to_fasta(
                self.opened_most_gff.name, self.args.reference,
                str(self.primary_fasta_path)))

    def get_link_region(self) -> None:
        """Get linked region for each primary region."""
        self.primary_sam_path = (self.tmp_floder_path /
                                 'primary.minimap2_reference.sam')
        self.output.info(f'Running minimap2 with {self.args.threads} threads')
        if self.args.verbose:
            self.output.info(splice_mapper(
                self.args.reference, str(self.primary_fasta_path),
                str(self.primary_sam_path), self.args.threads))
        else:
            self.opened_log_file.write('minimap2 info:\n')
            self.opened_log_file.write(splice_mapper(
                self.args.reference, str(self.primary_fasta_path),
                str(self.primary_sam_path), self.args.threads))
        self.primary_gff_path = (
            self.tmp_floder_path / 'primary.minimap2_reference.gff3')
        self.opened_primary_gff = self.primary_gff_path.open('w')
        self.output.info('Running cDNA_cupcake sam_to_gff3')
        if not self.args.verbose:
            sys.stdout = self.opened_log_file
            self.opened_log_file.write('cDNA_cupcake sam_to_gff3 info:\n')
        sam_to_gff(str(self.primary_sam_path),
                   self.opened_primary_gff, str(self.primary_fasta_path))
        if not self.args.verbose:
            sys.stdout = sys.__stdout__
        self.opened_primary_gff.close()
        self.opened_primary_gff = self.primary_gff_path.open('r')
        self.gene_linked_region = read_gff.read_gff(self.opened_primary_gff)
        self.link_id_list, self.linked_region_list = [], []
        for link_id, gene_id in enumerate(self.gene_linked_region.keys()):
            full_link_id = 'sgp_region' + str(link_id)
            self.link_id_list.append(full_link_id)
            self.gene_linked_region[gene_id].update_info_id(full_link_id)
            self.linked_region_list.append(self.gene_linked_region[gene_id])

    def process_links(self) -> None:
        """Using multiply threads process each linked region."""
        in_pool_threads = int(self.args.threads / len(self.link_id_list))
        in_pool_threads = in_pool_threads if in_pool_threads else 1
        out_pool_threads = int(self.args.threads / in_pool_threads)
        self.output.info(
            f'Calling {out_pool_threads} pools to index each region')
        process_link_args = []
        for link_id, linked_region in zip(
                self.link_id_list, self.linked_region_list):
            process_link_args.append((
                link_id, linked_region, self.tmp_floder_path,
                self.args.reference, self.args.input,
                self.args.fastx, in_pool_threads))
        with Pool(processes=out_pool_threads) as pool:
            pool.map(process_each_link, process_link_args)

    def process(self) -> None:
        """Call the index object."""
        self.output.info('Starting index Process')
        logger.debug('Starting index Process')
        self.prepare()
        self.check_bed()
        self.get_multimapped_reads()
        self.get_primary_region()
        self.collapse_primary_sam()
        self.get_most_supported_fasta()
        self.get_link_region()
        self.process_links()
        self.output.info('Completed index Process')
        logger.debug('Completed index Process')
