import gzip
import lzma
import bz2
from collections import defaultdict
from pathlib import Path
from typing import Union, Optional, Iterable, Tuple, Generator


openers = dict(
	gz=gzip.open,
	xz=lzma.open,
	bz2=bz2.open,
)
openers = defaultdict(lambda: open, **openers)


def transparent_open(
	file: Union[Path, str, Iterable[bytes]],
	mode: str = 'rt',
	*,
	ensure_parentdir: bool = False,
	encoding: Optional[str] = None,
	errors: Optional[str] = None,
	newline: Optional[str] = None,
	suffix: str = None
) -> Iterable[Union[str, bytes]]:
	"""
	Open potentially compressed file
	:param file: File path or file-like object opened in binary mode
	:param mode: Mode of the returned file
	:param ensure_parentdir: Ensure that the parent directory exists?
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
		if ensure_parentdir:
			path.parent.mkdir(parents=True, exist_ok=True)
	
	opener = openers[suffix]
	
	if not isinstance(file, str) and opener is open:
		return file
	else:
		try:
			return opener(file, mode, encoding=encoding, errors=errors, newline=newline)
		except TypeError as e:
			raise TypeError(f'Error in opener {opener}') from e


def parse_fq(line_header: str, line_seq: str, line_plus: str, line_qual: str) -> Tuple[str, str, str]:
	header = line_header.strip()
	assert header.startswith('@')
	seq_read = line_seq.strip()
	assert line_plus.startswith('+')
	seq_qual = line_qual.strip()
	return header, seq_read, seq_qual


def iter_fq(lines: Iterable[str]) -> Generator[Tuple[str, str, str], None, None]:
	line_it = iter(lines)
	for header in line_it:
		try:
			yield parse_fq(header, next(line_it), next(line_it), next(line_it))
		except StopIteration:
			raise IOError(f'Fastq file doesn’t contain new line after header {header}')


def read_fasta(filename: Union[Path, str]) -> Generator[Tuple[str, str], None, None]:
	with Path(filename).open() as f_bc:
		for header in f_bc:
			header = header.lstrip('>').strip()
			bc = next(f_bc).strip()
			yield header, bc
