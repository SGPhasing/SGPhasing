# -*- coding: utf-8 -*-
# Copyright 2021 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""SGPhasing is a Python3 package for ***.

This is a rough draft of the author's analysis of ***.
We will continue to add applications until we felt ready to say OK, and
we will release an official version through pip and bioconda.

For more in-depth instructions, see the SGPhasing manual, linked to from:

    https://github.com/SGPhasing/SGPhasing/wiki

If all else fails, feel free to post on GitHub Issue:

    https://github.com/SGPhasing/SGPhasing/issues

or contact Shang Xie and ask for help:

    xieshang0608@gmail.com
"""

from sys import exit, version_info

from SGPhasing import fullhelp_argumentparser
from SGPhasing.sys_output import Output

# version control
if version_info[0] == 3 and version_info[1] in (6, 7, 8, 9):
    pass
else:
    output = Output()
    output.error('Please run this script with Python version '
                 '3.6, 3.7, 3.8 or 3.9 and try again.')
    exit()


def main() -> None:
    """Creat subcommands and execute."""
    PARSER = fullhelp_argumentparser.FullHelpArgumentParser()
    SUBPARSER = PARSER.add_subparsers()
    INDEX = fullhelp_argumentparser.IndexArgs(
        SUBPARSER,
        'index',
        """Indexing with high quality full-length
            transcriptome sequences allele.""")
    PHASE = fullhelp_argumentparser.PhaseArgs(
        SUBPARSER,
        'phase',
        """.""")

    def bad_args(args) -> None:
        """Print help on bad arguments."""
        PARSER.print_help()
        exit()

    PARSER.set_defaults(func=bad_args)
    ARGUMENTS = PARSER.parse_args()
    ARGUMENTS.func(ARGUMENTS)


if __name__ == '__main__':
    main()
