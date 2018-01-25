import gzip
import lzma
import bz2
from collections import defaultdict
from pathlib import Path
from typing import Union, Optional, Iterable, Tuple, Generator, Pattern

from bartseq import BASES

openers = dict(
	gz=gzip.open,
	xz=lzma.open,
	bz2=bz2.open,
)
openers = defaultdict(lambda: open, **openers)


def transparent_open(
	file: Union[Path, str, Iterable[bytes]],
	mode: str='rt',
	*,
	encoding: Optional[str]=None,
	errors:   Optional[str]=None,
	newline:  Optional[str]=None,
	suffix: str=None
) -> Iterable[Union[str, bytes]]:
	"""
	Open potentially compressed file
	:param file: File path or file-like object opened in binary mode
	:param mode: Mode of the returned file
	:param encoding: Encoding of the text data the file decompresses to (or contains if the file is uncompressed)
	:param errors: See ``open``
	:param newline: See ``open``
	:param suffix: See ``open``
	:return: File-like object with decompressed data
	"""
	if isinstance(file, (str, Path)):
		path = Path(file)
		suffix = path.suffix[1:] if suffix is None else suffix
		file = str(file)
	
	opener = openers[suffix]
	
	if not isinstance(file, str) and opener is open:
		return file
	else:
		try:
			return opener(file, mode, encoding=encoding, errors=errors, newline=newline)
		except TypeError as e:
			raise TypeError(f'Error in opener {opener}') from e


def parse_fq(line_header: str, line_seq: str, line_plus: str, line_qual: str) -> Tuple[str, str, str]:
	header = line_header.rstrip('\n')
	assert header.startswith('@')
	seq_read = line_seq.rstrip('\n')
	assert line_plus.startswith('+')
	seq_qual = line_qual.rstrip('\n')
	return header, seq_read, seq_qual


def iter_fq(lines: Iterable[str]) -> Generator[Tuple[str, str, str], None, None]:
	line_it = iter(lines)
	while True:
		yield parse_fq(next(line_it), next(line_it), next(line_it), next(line_it))


def read_bcs(filename: str) -> Generator[Tuple[str, str], None, None]:
	with open(filename) as f_bc:
		header = next(f_bc)
		assert not all(c in BASES for field in header.split(' ') for c in field)
		for l in f_bc:
			id_, bc = l.rstrip('\n').split('\t')
			yield id_, bc
