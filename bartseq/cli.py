from typing import Dict

from .cli_helpers import CLI, DelegatingCLI
from .read_tagger.cli import ReadTaggerCLI


SUBCMDS = {
	'tag': ReadTaggerCLI(),
}  # type: Dict[str, CLI]


cli = DelegatingCLI(SUBCMDS)
