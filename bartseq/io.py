import gzip
import json
import lzma
import bz2
from collections import defaultdict
from pathlib import Path
from typing import Union, Optional, Iterable, Tuple, Generator

from . import BASES

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
	header = line_header.strip()
	assert header.startswith('@')
	seq_read = line_seq.strip()
	assert line_plus.startswith('+')
	seq_qual = line_qual.strip()
	return header, seq_read, seq_qual


def iter_fq(lines: Iterable[str]) -> Generator[Tuple[str, str, str], None, None]:
	line_it = iter(lines)
	while True:
		yield parse_fq(next(line_it), next(line_it), next(line_it), next(line_it))


def read_bcs(filename: Union[Path, str]) -> Generator[Tuple[str, str], None, None]:
	with Path(filename).open() as f_bc:
		for header in f_bc:
			header = header.lstrip('>').strip()
			bc = next(f_bc).strip()
			yield header, bc


def write_bc_table(path_bc_file: Union[Path, str], path_bc_table: Union[Path, str]):
	"""This is mainly independent of the rest, so do simple duplicate work to be able to create this separately"""
	from .read_tagger import get_tagger
	bc_table = get_tagger(read_bcs(path_bc_file)).get_barcode_table()
	with transparent_open(path_bc_table, 'wt') as f_bc:
		f_bc.write(bc_table)


def write_stats(stats_file: Union[Path, str], n_reads: int, n_both_regular: Optional[int], stats1: dict, stats2: Optional[dict]=None):
	with transparent_open(stats_file, 'wt') as f_s:
		stats = dict(n_reads=n_reads, read1=stats1)
		if n_both_regular is not None:
			stats['n_both_regular'] = n_both_regular
		if stats2:
			stats['read2'] = stats2
		json.dump(stats, f_s, indent='\t')
