# -*- coding: utf-8 -*-
# Copyright 2020 Shang Xie.
# All rights reserved.
#
# This file is part of the SGPhasing distribution and
# governed by your choice of the "SGPhasing License Agreement"
# or the "GNU General Public License v3.0".
# Please see the LICENSE file that should
# have been included as part of this package.
"""Represent an output.

What's here:

Format and display output.
--------------------------

Classes:
  - Output
"""

from platform import system

from rich import print as rprint


class Output():
    """Format and display output.

    Attributes:
        term_support_color: Term support color.
    """

    def __init__(self) -> None:
        """Initialize Output."""
        self.term_support_color = system() in ('Linux', 'Darwin')

    @staticmethod
    def __indent_text_block(text: str) -> str:
        """Indent a text block."""
        lines = text.splitlines()
        if len(lines) > 1:
            out = lines[0] + '\r\n'
            for i in range(1, len(lines) - 1):
                out = out + '        ' + lines[i] + '\r\n'
            out = out + '        ' + lines[-1]
            return out
        return text

    def info(self, text: str) -> None:
        """Format INFO text."""
        if text:
            trm = 'INFO    '
            if self.term_support_color:
                trm = (':information_source: ' +
                       '[bright_green]INFO[/bright_green]    ')
            rprint(trm + self.__indent_text_block(text))

    def warning(self, text: str) -> None:
        """Format WARNING text."""
        if text:
            trm = 'WARNING '
            if self.term_support_color:
                trm = ':warning: [bright_yellow]WARNING[/bright_yellow] '
            rprint(trm + self.__indent_text_block(text))

    def error(self, text: str) -> None:
        """Format ERROR text."""
        if text:
            trm = 'ERROR   '
            if self.term_support_color:
                trm = ':no_entry: [bright_red]ERROR[/bright_red]   '
            rprint(trm + self.__indent_text_block(text))
