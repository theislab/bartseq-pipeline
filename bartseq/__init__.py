from argparse import ArgumentParser
from pathlib import Path
from typing import Sequence, Callable, Generator, Tuple, Union, Iterable

import sys

import re

from .read_tagger import ReadTagger
from .io import transparent_open, openers


parser = ArgumentParser()
parser.add_argument(
	'in_file', nargs='?', default='-',
	help='File to read. supported compression: see --in-compression')
parser.add_argument(
	'--bc-file', '-b', required=True,
	help='Barcode file in the format ``<ID> <Sequence>`` (with header)')
parser.add_argument(
	'--in-compression', '-c', choices=openers.keys(),
	help='Specify compression if reading from stdin or a file with unusual compression')


def main(argv: Sequence[str]=None):
	if argv is None:
		argv = sys.argv[1:]
	
	args = parser.parse_args(argv)
	
	if args.in_file == '-':
		args.in_file = sys.stdin.buffer
	
	barcodes = {id: bc for id, bc in read_bcs(args.bc_file, get_pred(args.in_file))}
	print(barcodes)
	
	with transparent_open(args.in_file, suffix=args.in_compression) as f_in:
		tagger = ReadTagger(barcodes.values())
		for l in f_in:
			header = l.rstrip('\n')
			assert header.startswith('@')
			read = tagger.tag_read(next(f_in).rstrip('\n'))
			assert next(f_in).startswith('+')
			qual = next(f_in).rstrip('\n')
			
			lj = len(read.junk or [])
			ljb = lj + len(read.barcode or [])
			
			print(header)
			print(read)
			print('+')
			print(qual[:lj] or None, qual[lj:ljb] or None, qual[ljb:] or None)


def read_bcs(filename: str, pred: Callable[[str], bool]) -> Generator[Tuple[str, str], None, None]:
	with open(filename) as f_bc:
		header = next(f_bc)
		assert not all(c in 'ATGC' for c in header)
		for l in f_bc:
			id_, bc = l.rstrip('\n').split('\t')
			if pred(id_):
				yield id_, bc


def get_pred(filename: Union[Path, str, Iterable[str]]) -> Callable[[str], bool]:
	def pred(_):
		return True
	
	if not isinstance(filename, (Path, str)):
		return pred
	
	read = next(iter(re.finditer(r'_R(\d)_', str(filename))), None)
	if read:
		side = 'L' if read.group(1) == '1' else 'R'
		
		def pred(id_):
			return id_.startswith(side)
	
	return pred
