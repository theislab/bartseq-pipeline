import json
import sys
import re
from argparse import ArgumentParser, ArgumentError, Namespace
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Generator, Sequence, Union, Tuple, Iterable, Pattern, Optional

from tqdm import tqdm

from .logging import log, init_logging
from .read_tagger import ReadTagger, BASES
from .io import transparent_open, openers, iter_fq, read_bcs

parser = ArgumentParser()
parser.add_argument(
	'in_1', nargs='?', default='-',
	help='Read1 file to read from. Supported compression: see --in-compression')
parser.add_argument(
	'out_1', nargs='?', default='-',
	help='Read1 file to write to. Supported compression: see --out-compression')
arg_in_2 = parser.add_argument(
	'--in-2', nargs='?',
	help='Read2 file to read from. Supported compression: see --in-compression')
arg_out_2 = parser.add_argument(
	'--out-2', nargs='?',
	help='Read2 file to write to. Supported compression: see --out-compression')
parser.add_argument(
	'--stats-file', '-s', nargs='?', default='-',
	help='File to write final stats to (in JSON format)')
parser.add_argument(
	'--bc-file', '-b', required=True,
	help='Barcode file in the format ``<ID> <Sequence>`` (with header)')
parser.add_argument(
	'--bc-id-pat1', type=re.compile, default=re.compile(r'L\d+'),
	help='Read 1 Barcode ID regex. Default: ``L\d+`` for left barcodes')
parser.add_argument(
	'--bc-id-pat2', type=re.compile, default=re.compile(r'R\d+'),
	help='Read 2 Barcode ID regex. Default: ``R\d+`` for left barcodes')
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


def main(argv: Sequence[str]=None):
	init_logging()
	args = parse_args(argv)
	has_two_reads = bool(args.in_2)
	bcs_all = list(read_bcs(args.bc_file))
	
	tagger1 = get_tagger(1, bcs_all, args.bc_id_pat1, args.len_linker, args.len_primer)
	tagger2 = get_tagger(2, bcs_all, args.bc_id_pat2, args.len_linker, args.len_primer) if has_two_reads else None
	
	with tqdm(total=args.total) if args.total != 0 else ctx_dummy() as pb, \
		transparent_open(args.in_1,  'rt', suffix=args.in_compression)  as f_in_1, \
		transparent_open(args.out_1, 'wt', suffix=args.out_compression) as f_out_1, \
		transparent_open(args.in_2,  'rt', suffix=args.in_compression)  if has_two_reads else ctx_dummy() as f_in_2, \
		transparent_open(args.out_2, 'wt', suffix=args.out_compression) if has_two_reads else ctx_dummy() as f_out_2:
		
		n_reads = 0
		if has_two_reads:
			for r, (parts1, parts2) in enumerate(zip(iter_fq(f_in_1), iter_fq(f_in_2))):
				read1 = tagger1.tag_read(*parts1)
				read2 = tagger2.tag_read(*parts2)
				
				try:
					if read1.barcode and read2.barcode:
						f_out_1.write(str(read1))
						f_out_2.write(str(read2))
				except BrokenPipeError:
					break
				
				update_pb(pb, tagger1, r)
				n_reads = r
				# TODO: maybe both?
		else:
			for r, (header, seq_read, seq_qual) in enumerate(iter_fq(f_in_1)):
				read = tagger1.tag_read(header, seq_read, seq_qual)
				
				try:
					f_out_1.write(str(read))
				except BrokenPipeError:
					break
				
				update_pb(pb, tagger1, r)
				n_reads = r
		
		if pb: pb.close()
	
	with transparent_open(args.stats_file, 'wt') as f_s:
		stats = dict(n_reads=n_reads, read1=tagger1.stats)
		if has_two_reads:
			stats['read2'] = tagger2.stats
		json.dump(stats, f_s, indent='\t')


def parse_args(argv: Sequence[str]=None) -> Union[Namespace, Any]:
	if argv is None:
		argv = sys.argv[1:]
	args = parser.parse_args(argv)
	
	if args.in_1 == '-': args.in_1 = sys.stdin
	if args.in_2 == '-': args.in_2 = sys.stdin
	
	if args.out_1 == '-': args.out_1 = sys.stdout
	else: Path(args.out_1).parent.mkdir(parents=True, exist_ok=True)
	if args.out_2 == '-': args.out_2 = sys.stdout
	else: Path(args.out_2).parent.mkdir(parents=True, exist_ok=True)
	
	if args.in_1 == args.in_2:
		raise ArgumentError(arg_in_2, 'Cannot parse both reads from the same file.')
	if args.out_1 == args.out_2:
		raise ArgumentError(arg_out_2, 'Cannot write both reads from the same file.')
	
	if bool(args.in_2) != bool(args.out_2):
		raise ArgumentError(arg_in_2, 'You need to specify both or none of --in-2 and --out-2.')
	
	return args


def get_tagger(r: int, bcs_all: Iterable[Tuple[str, str]], pat: Pattern, len_linker: int, len_primer: int):
	bc_to_id = {bc: id_ for id_, bc in bcs_all if pat.match(id_)}
	log.info(f'Using read {r} barcodes: {bc_to_id}')
	return ReadTagger(bc_to_id, len_linker, len_primer)


def update_pb(pb: Optional[tqdm], tgr: ReadTagger, i: int):
	if not pb: return
	pb.update(1)
	if i % 10000 == 0:
		pb.set_description(', '.join(f'{k}: {v}' for k, v in tgr.stats.items()), refresh=False)


@contextmanager
def ctx_dummy() -> Generator[None, None, None]:
	yield None
