import gzip
import lzma
import bz2
from collections import defaultdict
from pathlib import Path
from typing import Union, Optional

openers = dict(
	gz=gzip.open,
	xz=lzma.open,
	bz2=bz2.open,
)
openers = defaultdict(lambda: open, **openers)


def transparent_open(
	filename: Union[Path, str],
	mode: str='rt',
	*,
	encoding: Optional[str]=None,
	errors:   Optional[str]=None,
	newline:  Optional[str]=None,
	suffix: str=None
):
	filename = Path(filename)
	
	suffix = filename.suffix[1:] if suffix is None else suffix
	return openers[suffix](str(filename), mode, encoding=encoding, errors=errors, newline=newline)
