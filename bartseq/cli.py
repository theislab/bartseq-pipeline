from typing import Dict

from .cli_helpers import CLI, DelegatingCLI
from .read_tagger.cli import ReadTaggerCLI
from .fastq_browser.cli import FastqBrowserCLI


SUBCMDS = {
	'tag': ReadTaggerCLI(),
	'browse': FastqBrowserCLI(),
}  # type: Dict[str, CLI]


cli = DelegatingCLI(SUBCMDS)
