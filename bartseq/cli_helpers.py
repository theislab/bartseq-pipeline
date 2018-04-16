from argparse import ArgumentParser, Namespace, _SubParsersAction
from typing import Sequence

import sys


class AbstractAttribute:
	def __get__(self, obj, type):
		raise NotImplemented


class CLI:
	COMMAND = AbstractAttribute()
	
	@staticmethod
	def populate_parser(parser: ArgumentParser) -> ArgumentParser:
		return parser or ArgumentParser()
	
	@staticmethod
	def modify_args(parser: ArgumentParser, args: Namespace) -> Namespace:
		return args
	
	@classmethod
	def register_as_subcmd(cls, subcommands: _SubParsersAction):
		subparser = cls.populate_parser(subcommands.add_parser(cls.COMMAND))
		
		def modify_and_run(args: Namespace):
			return cls.run(cls.modify_args(subparser, args))
		
		subparser.set_defaults(func=modify_and_run)
		return subparser
	
	@classmethod
	def run_as_main(cls, argv: Sequence[str] = None):
		if argv is None:
			argv = sys.argv[1:]
		parser = ArgumentParser()
		args = parser.parse_args(argv)
		args = cls.modify_args(parser, args)
		cls.run(args)
	
	@staticmethod
	def run(args: Namespace):
		raise NotImplemented
