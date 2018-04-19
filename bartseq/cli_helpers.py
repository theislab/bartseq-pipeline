from argparse import ArgumentParser, Namespace, _SubParsersAction
from typing import Sequence, Dict

import sys


class AbstractAttribute:
	def __get__(self, obj, type):
		raise NotImplemented


class CLI:
	@staticmethod
	def populate_parser(parser: ArgumentParser) -> ArgumentParser:
		return parser or ArgumentParser()
	
	@staticmethod
	def modify_args(parser: ArgumentParser, args: Namespace) -> Namespace:
		return args
	
	def run_as_main(self, argv: Sequence[str] = None):
		if argv is None:
			argv = sys.argv[1:]
		parser = self.populate_parser(ArgumentParser())
		args = parser.parse_args(argv)
		args = self.modify_args(parser, args)
		self.run(parser, args)
	
	@staticmethod
	def run(parser: ArgumentParser, args: Namespace):
		raise NotImplemented


class DelegatingCLI(CLI):
	def __init__(self, subcmds: Dict[str, CLI]):
		self.subcmds = subcmds
	
	@staticmethod
	def _register_subcommand(subcmd: CLI, parser: ArgumentParser) -> ArgumentParser:
		subparser = subcmd.populate_parser(parser)
		
		def modify_and_run(args: Namespace):
			return subcmd.run(subparser, subcmd.modify_args(subparser, args))
		
		subparser.set_defaults(func=modify_and_run)
		return subparser
	
	def populate_parser(self, parser: ArgumentParser) -> ArgumentParser:
		subparsers = parser.add_subparsers()  # type: _SubParsersAction
		
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
