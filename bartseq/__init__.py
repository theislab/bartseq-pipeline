from argparse import ArgumentParser
from typing import Sequence

import sys

from .read_tagger import ReadTagger
from .io import transparent_open, openers


parser = ArgumentParser()
parser.add_argument('in_file', default='-', help='File to read. supported compression: see --in-compression')
parser.add_argument('--in-compression', '-i', choices=openers.keys())


def main(argv: Sequence[str]=None):
	if argv is None:
		argv = sys.argv[1:]
	
	args = parser.parse_args(argv)
	if args.in_file == '-':
		raise NotImplemented
	
	with transparent_open(args.in_file, suffix=args.in_compression) as f_in:
		tagger = ReadTagger(['TODO'])  # TODO
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
