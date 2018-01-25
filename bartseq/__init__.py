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
from .io import transparent_open, openers, iter_fq
from .logging import log, init_logging

parser = ArgumentParser()
parser.add_argument(
	'in_file', nargs='?', default='-',
	help='File to read. Supported compression: see --in-compression')
parser.add_argument(
	'out_file', nargs='?', default='-',
	help='File to write to. Supported compression: see --out-compression')
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
def ctx_dummy() -> Generator[None, None, None]:
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
	else:
		Path(args.out_file).parent.mkdir(parents=True, exist_ok=True)
	
	id_to_bc = {id_: bc for id_, bc in read_bcs(args.bc_file, args.bc_id_pat)}
	bc_to_id = {bc: id_ for id_, bc in id_to_bc.items()}
	log.info(f'Using barcodes: {id_to_bc}')
	tagger = ReadTagger(bc_to_id.keys(), args.len_linker, args.len_primer)
	
	with tqdm(total=args.total*4) if args.total != 0 else ctx_dummy() as pb, \
		transparent_open(args.in_file,  'rt', suffix=args.in_compression) as f_in, \
		transparent_open(args.out_file, 'wt', suffix=args.out_compression) as f_out:
		
		for l, (header, seq_read, seq_qual) in enumerate(iter_fq(f_in)):
			read = tagger.tag_read(header, seq_read, seq_qual)
			
			try:
				f_out.write(str(read))
			except BrokenPipeError:
				break
			
			if pb:
				pb.update(1)
				if l % 10000 == 0:
					pb.set_description(', '.join(f'{k}: {v}' for k, v in tagger.stats.items()), refresh=False)
		
		if pb: pb.close()
	
	with transparent_open(args.stats_file, 'wt') as f_s:
		json.dump(tagger.stats, f_s)


def read_bcs(filename: str, id_pat: Pattern) -> Generator[Tuple[str, str], None, None]:
	with open(filename) as f_bc:
		header = next(f_bc)
		assert not all(c in BASES for field in header.split(' ') for c in field)
		for l in f_bc:
			id_, bc = l.rstrip('\n').split('\t')
			if id_pat.match(id_):
				yield id_, bc
