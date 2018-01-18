import json
from argparse import ArgumentParser
from contextlib import contextmanager, ExitStack
from pathlib import Path
from textwrap import dedent
from typing import Sequence, Generator, Tuple, Pattern

import sys

import re

from tqdm import tqdm

from .read_tagger import ReadTagger, BASES
from .io import transparent_open, openers
from .logging import log, init_logging

parser = ArgumentParser()
parser.add_argument(
	'in_file', nargs='?', default='-',
	help='File to read. Supported compression: see --in-compression')
parser.add_argument(
	'out_file', nargs='?', default='-',
	help='File to write to. If it contains “{}”, demultiplexing is used. Supported compression: see --out-compression')
parser.add_argument(
	'--stats-file', '-s', nargs='?', default='-',
	help='File to write final stats to (in JSON format)')
parser.add_argument(
	'--bc-file', '-b', required=True,
	help='Barcode file in the format ``<ID> <Sequence>`` (with header)')
parser.add_argument(
	'--bc-id-pat', '-r', type=re.compile, required=True,
	help='Barcode ID regex. e.g. ``L\d+`` for left barcodes')
parser.add_argument(
	'--total', '-t', type=int, default=0,
	help='Number of fastq records in file. “0” means no progressbar')
parser.add_argument(
	'--len-primer', '-p', type=int, default=27,
	help='Primer length for stats')
parser.add_argument(
	'--len-linker', '-l', type=int, default=10,
	help='Linker length to cut out')
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
	init_logging()
	if argv is None:
		argv = sys.argv[1:]
	
	args = parser.parse_args(argv)
	
	if args.in_file == '-':
		args.in_file = sys.stdin
	if args.out_file == '-':
		args.out_file = sys.stdout
		multiple_outputs = False
	else:
		Path(args.out_file).parent.mkdir(parents=True, exist_ok=True)
		multiple_outputs = '{}' in args.out_file
	
	barcodes = {id_: bc for id_, bc in read_bcs(args.bc_file, args.bc_id_pat)}
	log.info(f'Using barcodes: {barcodes}')
	tagger = ReadTagger(barcodes.values(), args.len_linker)
	
	stats = dict(
		n_only_primer=0,
		n_multiple_bcs=0,
		n_no_barcode=0,
		n_regular=0,
		n_junk=0,
	)
	
	with tqdm(total=args.total*4) if args.total != 0 else ctx_dummy() as pb, \
		transparent_open(args.in_file,  'rt', suffix=args.in_compression) as f_in, \
		ExitStack() as stack:
		
		fs_out = {}
		
		if multiple_outputs:
			fs_out[None] = transparent_open(args.out_file.format('unmapped'), 'wt', suffix=args.out_compression)
			for bc in barcodes.values():
				fs_out[bc] = transparent_open(args.out_file.format(bc), 'wt', suffix=args.out_compression)
		else:
			fs_out[None] = transparent_open(args.out_file, 'wt', suffix=args.out_compression)
		
		for l, line in enumerate(f_in):
			header = line.rstrip('\n')
			assert header.startswith('@')
			read = tagger.tag_read(next(f_in).rstrip('\n'))
			assert next(f_in).startswith('+')
			qual = read.cut_seq(next(f_in))
			
			is_primer = read.is_just_primer(args.len_primer)
			
			if is_primer:                      stats['n_only_primer'] += 1
			if read.has_multiple_barcodes:     stats['n_multiple_bcs'] += 1
			if not read.barcode:               stats['n_no_barcode'] += 1
			if read.junk:                      stats['n_junk'] += 1
			if read.barcode and not is_primer: stats['n_regular'] += 1
			
			tagline = ' '.join([
				f'barcode={read.barcode}',
				f'linker={read.linker}',
				f'multi-bc={read.has_multiple_barcodes}',
				f'just-primer={is_primer}',
				f'other-bcs={",".join(read.other_barcodes) or None}',
				f'junk={read.junk}',
			])
			
			f_out = fs_out[read.barcode if multiple_outputs else None]
			try:
				f_out.write(dedent(f'''\
					{header} {tagline}
					{read.amplicon}
					+
					{qual}
				'''))
				if pb:
					pb.update(1)
					if l % 10000 == 0:
						pb.set_description(', '.join(f'{k}: {v}' for k, v in stats.items()))
			except BrokenPipeError:
				break
		
		if pb: pb.close()
	
	with transparent_open(args.stats_file, 'wt') as f_s:
		json.dump(stats, f_s)


def read_bcs(filename: str, id_pat: Pattern) -> Generator[Tuple[str, str], None, None]:
	with open(filename) as f_bc:
		header = next(f_bc)
		assert not all(c in BASES for field in header.split(' ') for c in field)
		for l in f_bc:
			id_, bc = l.rstrip('\n').split('\t')
			if id_pat.match(id_):
				yield id_, bc
