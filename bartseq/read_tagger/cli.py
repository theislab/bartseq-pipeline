from argparse import ArgumentParser, Action, Namespace, ArgumentError
from pathlib import Path

import sys

from . import defaults
from ..io import openers
from ..cli_helpers import CLI


class ReadTaggerCLI(CLI):
	COMMAND = 'tag'
	
	@staticmethod
	def populate_parser(parser: ArgumentParser) -> ArgumentParser:
		parser.add_argument(
			'in_1', nargs='?', default='-',
			help='Read1 file to read from. Supported compression: see --in-compression')
		parser.add_argument(
			'out_1', nargs='?', default='-',
			help='Read1 file to write to. Supported compression: see --out-compression')
		parser.add_argument(
			'--in-2', nargs='?',
			help='Read2 file to read from. Supported compression: see --in-compression')
		parser.add_argument(
			'--out-2', nargs='?',
			help='Read2 file to write to. Supported compression: see --out-compression')
		parser.add_argument(
			'--bc-file', '-b', required=True,
			help='Barcode file in the format ``<ID> <Sequence>`` (with header)')
		parser.add_argument(
			'--stats-file', '-s', required=True,
			help='File to write final stats to (in JSON format)')
		parser.add_argument(
			'--bc-table', '-B', nargs='?',
			help='File name for the HTML table of barcode mismatches')
		parser.add_argument(
			'--total', '-t', type=int, default=defaults.total,
			help='Number of fastq records in file. “0” means no progressbar')
		parser.add_argument(
			'--len-primer', '-p', type=int, default=defaults.len_primer,
			help='Primer length for stats')
		parser.add_argument(
			'--len-linker', '-l', type=int, default=defaults.len_linker,
			help='Linker length to cut out')
		parser.add_argument(
			'--in-compression', '-i', choices=openers.keys(),
			help='Specify compression if reading from stdin or a file with unusual suffix')
		parser.add_argument(
			'--out-compression', '-o', choices=openers.keys(),
			help='Specify compression if writing to stdout or a file with unusual suffix')
		
		return parser
	
	@staticmethod
	def modify_args(parser: ArgumentParser, args: Namespace) -> Namespace:
		if args.in_1 == '-': args.in_1 = sys.stdin
		if args.in_2 == '-': args.in_2 = sys.stdin
		
		if args.out_1 == '-': args.out_1 = sys.stdout
		else: Path(args.out_1).parent.mkdir(parents=True, exist_ok=True)
		if args.out_2 == '-': args.out_2 = sys.stdout
		else: Path(args.out_2).parent.mkdir(parents=True, exist_ok=True)
		
		def find_action(dest: str) -> Action:
			return next(action for action in parser._actions if action.dest == dest)
		
		if args.in_1 == args.in_2:
			raise ArgumentError(find_action('in_2'), 'Cannot parse both reads from the same file.')
		if args.out_1 == args.out_2:
			raise ArgumentError(find_action('out_2'), 'Cannot write both reads from the same file.')
		
		if bool(args.in_2) != bool(args.out_2):
			raise ArgumentError(find_action('in_2'), 'You need to specify both or none of --in-2 and --out-2.')
		
		return args
	
	@staticmethod
	def run(args: Namespace):
		from .main import run
		run(**vars(args))
