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
from subprocess import Popen

from SGPhasing.reader.read_xam import check_index
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
        minimap2_index = Path(self.args.ref + '.mmi')
        if not minimap2_index.exists():
            with Popen(['minimap2', '-k', '17', '-I', '18G', '-x', 'splice',
                        '-d', f'{self.args.ref}.mmi', self.args.ref]):
                self.output.info('preparing genome index for minimap2.')
        check_index(self.args.input)

    def process(self) -> None:
        """Call the index object."""
        logger.debug('Starting index Process')
        self.prepare()
        logger.debug('Completed index Process')
