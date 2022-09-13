"""General utilities to parse both input/out files."""

__all__ = ['any_char', 'integer', 'natural', 'parse_file', 'parse_section',
           'skipany_char', 'skip_line', 'skip_supress', 'try_search_pattern']


import re
from pathlib import Path
from typing import Optional as Optional_
from typing import Union


from pyparsing import (Combine, Literal, Optional,
                       ParseException, ParserElement, ParseResults, Regex,
                       SkipTo, Suppress, Word, nums)

PathLike = Union[Path, str]

# Literals
minus_or_plus = Literal('+') | Literal('-')

# Parsing Floats
natural = Word(nums)
integer = Combine(Optional(minus_or_plus) + natural)
float_number = Regex(r'(\-)?\d+(\.)(\d*)?([eE][\-\+]\d+)?')

float_number_dot = Regex(r'(\-)?(\d+)?(\.)(\d*)?([eE][\-\+]\d+)?')

# Parse Utilities


def skip_supress(z: str) -> ParserElement:
    """Skip until `z` and suppress the skipped values."""
    return Suppress(SkipTo(z))


any_char = Regex('.')
skipany_char = Suppress(any_char)
skip_line = Suppress(skip_supress('\n'))

# Generic Functions


def parse_file(p: ParserElement, file_name: PathLike) -> ParseResults:
    """Apply parser `p` on file `file_name`."""
    try:
        return p.parseFile(file_name)
    except ParseException as ex:
        raise ParseException(f"Error Trying to parse: {p} in file: {file_name}") from ex


def parse_section(start: str, end: str) -> ParserElement:
    """Read the lines from `start` to `end`."""
    s = Literal('{}'.format(start))
    e = Literal('{}'.format(end))

    return Suppress(SkipTo(s)) + skip_line + SkipTo(e)


def try_search_pattern(pat: str, file_name: PathLike) -> Optional_[str]:
    """Search for an specific pattern in  a file."""
    try:
        with open(file_name, 'r') as f:
            for line in f:
                if re.search(pat, line):
                    return line
            else:
                return None
    except FileNotFoundError:
        msg2 = f'There is not a file: {file_name}\n'
        raise RuntimeError(msg2)
