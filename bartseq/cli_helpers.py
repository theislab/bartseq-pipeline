import sys
from argparse import ArgumentParser, Namespace, _SubParsersAction
from pathlib import Path
from typing import Sequence, Dict, Union, TextIO


class AbstractAttribute:
	def __get__(self, obj, type):
		raise NotImplemented


def t_in_file(f: Union[Path, str, TextIO]) -> Union[Path, str, TextIO]:
	if f == '-': return sys.stdin
	return f


def t_out_file(f: Union[Path, str, TextIO]) -> Union[Path, str, TextIO]:
	if f == '-': return sys.stdout
	return f


class CLI:
	@staticmethod
	def populate_parser(parser: ArgumentParser) -> ArgumentParser:
		return parser or ArgumentParser()
	
	@staticmethod
	def check_args(parser: ArgumentParser, args: Namespace):
		"""Override to check and optionally modify arguments."""
		pass
	
	@staticmethod
	def run(parser: ArgumentParser, args: Namespace):
		raise NotImplemented
	
	def check_and_run(self, parser: ArgumentParser, args: Namespace):
		self.check_args(parser, args)
		self.run(parser, args)
	
	def run_as_main(self, argv: Sequence[str] = None):
		if argv is None:
			argv = sys.argv[1:]
		parser = self.populate_parser(ArgumentParser())
		args = parser.parse_args(argv)
		self.check_and_run(parser, args)


class DelegatingCLI(CLI):
	def __init__(self, subcmds: Dict[str, CLI]):
		self.subcmds = subcmds
	
	@staticmethod
	def _register_subcommand(subcmd: CLI, parser: ArgumentParser) -> ArgumentParser:
		subparser = subcmd.populate_parser(parser)
		
		def check_and_run(args: Namespace):
			subcmd.check_and_run(subparser, args)
		
		subparser.set_defaults(func=check_and_run)
		return subparser
	
	def populate_parser(self, parser: ArgumentParser) -> ArgumentParser:
		subparsers: _SubParsersAction = parser.add_subparsers()
		
		for name, subcmd in self.subcmds.items():
			self._register_subcommand(subcmd, subparsers.add_parser(name))
		
		return parser
	
	def run(self, parser: ArgumentParser, args: Namespace):
		if hasattr(args, 'func'):
			return args.func(args)
		elif not vars(args):
			parser.print_help()
		else:
			assert False, 'This should never happen'
