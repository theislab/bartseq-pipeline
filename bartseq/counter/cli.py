from argparse import ArgumentParser, Namespace
from pathlib import Path

from .main import main
from ..cli_helpers import CLI, clean_kbdinterrupt, suggest_library


class CounterCLI(CLI):
	@staticmethod
	def populate_parser(parser: ArgumentParser) -> ArgumentParser:
		parser.add_argument(
			'data_dir', type=Path,
			help='Data directory to read from. Needs to have the directories “./process/{3-tagged,4-mapped}” filled.')
		parser.add_argument(
			'library', nargs='?', default=None, help=(
				'Library name. E.g. “Lib1_S1_L001” for input files named “Lib1_S1_L001_R{12}_001.fastq.gz”. '
				'Omittable if only one library exists.'))
		parser.add_argument(
			'--no-mismatch', dest='allow_mismatch', default=True, action='store_false',
			help='Ignore barcodes with mismatches while counting.')
		parser.add_argument(
			'--both', default=None, action='store_true',
			help='Print the count results for both to stdout. Default: Write to “./process/5-counts” instead')
		parser.add_argument(
			'--one', default=None, action='store_true',
			help='Print the count results for one to stdout. Default: Write to “./process/5-counts” instead')
		return parser
	
	@staticmethod
	def check_args(parser: ArgumentParser, args: Namespace):
		if (
			(args.both and args.one is not None) or
			(args.one and args.both is not None)
		): parser.error('Cannot specify --both and --one')
		if args.one:
			args.both = not args.one
		del args.one
		args.library = suggest_library(args.data_dir, args.library, parser.error)
	
	@staticmethod
	@clean_kbdinterrupt
	def run(parser: ArgumentParser, args: Namespace):
		kwargs = vars(args)
		del kwargs['func']
		main(**kwargs)


cli = CounterCLI()
