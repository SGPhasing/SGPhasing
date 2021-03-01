# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing.processor gatk4.

Functions:
  - create_sequence_dictionary
  - add_or_replace_read_groups
  - left_align_indels
  - haplotype_caller
"""

from subprocess import PIPE, Popen


def create_sequence_dictionary(reference: str) -> str:
    """Call gatk4 for creating sequence dictionary.

    Args:
        reference (str): input reference fasta file path string.

    Returns:
        (str): gatk4 print.
    """
    with Popen(['gatk', 'CreateSequenceDictionary', '--REFERENCE', reference],
               stdout=PIPE) as proc:
        return proc.stdout.read().decode('utf-8')


def add_or_replace_read_groups(input_bam: str,
                               output_bam: str,
                               rglb: str,
                               rgpl: str,
                               rgpu: str,
                               rgsm: str) -> str:
    """Call gatk4 for addind or replacing read groups.

    Args:
        input_bam (str): input bam file path string.
        output_bam (str): output bam file path string.
        rglb (str): read-group library.
        rgpl (str): read-group platform (e.g. ILLUMINA, SOLID).
        rgpu (str): read-group platform unit (eg. run barcode).
        rgsm (str): read-group sample name.

    Returns:
        (str): gatk4 print.
    """
    with Popen(['gatk', 'AddOrReplaceReadGroups', '--INPUT', input_bam,
                '--OUTPUT', output_bam, '--RGLB', rglb, '--RGPL', rgpl,
                '--RGPU', rgpu, '--RGSM', rgsm], stdout=PIPE) as proc:
        return proc.stdout.read().decode('utf-8')


def left_align_indels(input_bam: str, output_bam: str, reference: str) -> str:
    """Call gatk4 for left aligning indels.

    Args:
        input_bam (str): input bam file path string.
        output_bam (str): output bam file path string.
        reference (str): input reference fasta file path string.

    Returns:
        (str): gatk4 print.
    """
    with Popen(['gatk', 'LeftAlignIndels',
                '--disable-tool-default-read-filters', 'true',
                '--input', input_bam, '--output', output_bam,
                '--reference', reference], stdout=PIPE) as proc:
        return proc.stdout.read().decode('utf-8')


def haplotype_caller(input_bam: str,
                     output_vcf: str,
                     reference: str,
                     max_reads: int,
                     min_quality: int,
                     ploidy: int,
                     threads: int = 1) -> str:
    """Call gatk4 for calling haplotype.

    Args:
        input_bam (str): input bam file path string.
        output_vcf (str): output vcf file path string.
        reference (str): input reference fasta file path string.
        max_reads (int): max reads number for --max-reads-per-alignment-start.
        min_quality (int): min quality number for --min-base-quality-score.
        ploidy (int): ploidy number for --sample-ploidy.
        threads (int): threads using for IntelPairHmm.

    Returns:
        (str): gatk4 print.
    """
    args = [
        'gatk', '--java-options', "'-DGATK_STACKTRACE_ON_USER_EXCEPTION=true'",
        'HaplotypeCaller', '--input', input_bam, '--output', output_vcf,
        '--reference', reference,
        '--max-reads-per-alignment-start', str(max_reads),
        '--min-base-quality-score', str(min_quality),
        '--native-pair-hmm-threads', str(threads),
        '--sample-ploidy', str(ploidy),
        '--do-not-run-physical-phasing', 'true',
        '--dont-use-soft-clipped-bases', 'true',
        '--max-alternate-alleles', str(ploidy*2),
        '--pcr-indel-model', 'AGGRESSIVE']
    with Popen(' '.join(args), stdout=PIPE, shell=True) as proc:
        return proc.stdout.read().decode('utf-8')
