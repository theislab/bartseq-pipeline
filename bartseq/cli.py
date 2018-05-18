from typing import Dict

from .cli_helpers import CLI, DelegatingCLI
from .read_tagger.cli import ReadTaggerCLI
from .fastq_browser.cli import FastqBrowserCLI


SUBCMDS: Dict[str, CLI] = {
	'tag': ReadTaggerCLI(),
	'browse': FastqBrowserCLI(),
}


cli = DelegatingCLI(SUBCMDS)
