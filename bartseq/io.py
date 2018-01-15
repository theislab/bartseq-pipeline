import gzip
import lzma
import bz2
from collections import defaultdict
from pathlib import Path
from typing import Union, Optional, Iterable

openers = dict(
	gz=gzip.open,
	xz=lzma.open,
	bz2=bz2.open,
)
openers = defaultdict(lambda: open, **openers)


def transparent_open(
	file: Union[Path, str, Iterable[bytes]],
	mode: str='rt',
	*,
	encoding: Optional[str]=None,
	errors:   Optional[str]=None,
	newline:  Optional[str]=None,
	suffix: str=None
) -> Iterable[Union[str, bytes]]:
	"""
	Open potentially compressed file
	:param file: File path or file-like object opened in binary mode
	:param mode: Mode of the returned file
	:param encoding: Encoding of the text data the file decompresses to (or contains if the file is uncompressed)
	:param errors: See ``open``
	:param newline: See ``open``
	:param suffix: See ``open``
	:return: File-like object with decompressed data
	"""
	if isinstance(file, (str, Path)):
		path = Path(file)
		suffix = path.suffix[1:] if suffix is None else suffix
		file = str(file)
	
	opener = openers[suffix]
	
	if not isinstance(file, str) and opener is open:
		return file
	else:
		try:
			return opener(file, mode, encoding=encoding, errors=errors, newline=newline)
		except TypeError as e:
			raise TypeError(f'Error in opener {opener}') from e
