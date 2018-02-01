from argparse import ArgumentParser, Action, Namespace, ArgumentError
from pathlib import Path
from typing import Sequence, Any, Union

import sys

from .io import openers


parser = ArgumentParser()
parser.add_argument(
	'in_1', nargs='?', default='-',
	help='Read1 file to read from. Supported compression: see --in-compression')
parser.add_argument(
	'out_1', nargs='?', default='-',
	help='Read1 file to write to. Supported compression: see --out-compression')
arg_in_2 = parser.add_argument(
	'--in-2', nargs='?',
	help='Read2 file to read from. Supported compression: see --in-compression')  # type: Action
arg_out_2 = parser.add_argument(
	'--out-2', nargs='?',
	help='Read2 file to write to. Supported compression: see --out-compression')  # type: Action
parser.add_argument(
	'--bc-file', '-b', required=True,
	help='Barcode file in the format ``<ID> <Sequence>`` (with header)')
parser.add_argument(
	'--stats-file', '-s', required=True,
	help='File to write final stats to (in JSON format)')
parser.add_argument(
	'--bc-table', '-B', nargs='?',
	help='File name for the HTML table of barcode mismatches')
arg_total = parser.add_argument(
	'--total', '-t', type=int, default=0,
	help='Number of fastq records in file. “0” means no progressbar')  # type: Action
arg_len_primer = parser.add_argument(
	'--len-primer', '-p', type=int, default=27,
	help='Primer length for stats')  # type: Action
arg_len_linker = parser.add_argument(
	'--len-linker', '-l', type=int, default=10,
	help='Linker length to cut out')  # type: Action
parser.add_argument(
	'--in-compression', '-i', choices=openers.keys(),
	help='Specify compression if reading from stdin or a file with unusual suffix')
parser.add_argument(
	'--out-compression', '-o', choices=openers.keys(),
	help='Specify compression if writing to stdout or a file with unusual suffix')


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


def run_cli(argv: Sequence[str]=None):
	from .main import run
	run(**vars(parse_args(argv)))
