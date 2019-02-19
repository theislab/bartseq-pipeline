from typing import Dict

from .cli_helpers import CLI, DelegatingCLI
from .read_tagger.cli import ReadTaggerCLI
from .fastq_browser.cli import FastqBrowserCLI
from .counter.cli import CounterCLI


SUBCMDS: Dict[str, CLI] = {
	'tag': ReadTaggerCLI(),
	'browse': FastqBrowserCLI(),
	'count': CounterCLI(),
}


cli = DelegatingCLI(SUBCMDS)
