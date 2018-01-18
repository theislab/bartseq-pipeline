from argparse import ArgumentParser
from contextlib import contextmanager
from pathlib import Path
from typing import Sequence, Callable, Generator, Tuple, Union, Iterable

import sys

import re

from tqdm import tqdm

from .read_tagger import ReadTagger, BASES
from .io import transparent_open, openers

parser = ArgumentParser()
parser.add_argument(
	'in_file', nargs='?', default='-',
	help='File to read. supported compression: see --in-compression')
parser.add_argument(
	'out_file', nargs='?', default='-',
	help='File to write to. supported compression: see --out-compression')
parser.add_argument(
	'--bc-file', '-b', required=True,
	help='Barcode file in the format ``<ID> <Sequence>`` (with header)')
parser.add_argument(
	'--total', '-t', type=int, default=0,
	help='Number of fastq records in file. “0” means no progressbar')
parser.add_argument(
	'--len-primer', '-p', type=int, default=22,
	help='Primer length for stats')
parser.add_argument(
	'--in-compression', '-i', choices=openers.keys(),
	help='Specify compression if reading from stdin or a file with unusual suffix')
parser.add_argument(
	'--out-compression', '-o', choices=openers.keys(),
	help='Specify compression if writing to stdout or a file with unusual suffix')


@contextmanager
def ctx_dummy():
	yield None


def main(argv: Sequence[str]=None):
	if argv is None:
		argv = sys.argv[1:]
	
	args = parser.parse_args(argv)
	
	if args.in_file == '-':
		args.in_file = sys.stdin
	if args.out_file == '-':
		args.out_file = sys.stdout
	
	barcodes = {id: bc for id, bc in read_bcs(args.bc_file, get_pred(args.in_file))}
	tagger = ReadTagger(barcodes.values())
	
	with tqdm(total=args.total) if args.total != 0 else ctx_dummy() as pb, \
		transparent_open(args.in_file,  'rt', suffix=args.in_compression) as f_in, \
		transparent_open(args.out_file, 'wt', suffix=args.out_compression) as f_out:
		
		for l in f_in:
			header = l.rstrip('\n')
			assert header.startswith('@')
			read = tagger.tag_read(next(f_in).rstrip('\n'))
			assert next(f_in).startswith('+')
			qual = next(f_in).rstrip('\n')
			
			ljb = len(read.junk or []) + len(read.barcode or [])
			ljba = ljb + len(read.amplicon or [])
			
			try:
				f_out.write(f'''\
{header} barcode={read.barcode} \
multi-bc={read.has_multiple_barcodes} \
just-primer={read.is_just_primer(args.len_primer)} \
other-bcs={",".join(read.other_barcodes) or None} \
junk={read.junk}
{read.amplicon}
+
{qual[ljb:ljba]}
''')
				if pb: pb.update(1)
			except BrokenPipeError:
				break
		
		if pb: pb.close()


def read_bcs(filename: str, pred: Callable[[str], bool]) -> Generator[Tuple[str, str], None, None]:
	with open(filename) as f_bc:
		header = next(f_bc)
		assert not all(c in BASES for field in header.split(' ') for c in field)
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
