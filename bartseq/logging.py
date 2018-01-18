from logging import (
	Logger, basicConfig,
	DEBUG, INFO, WARNING, ERROR, CRITICAL,
)

import sys

log = Logger(__name__)


def init_logging():
	basicConfig(level=INFO, stream=sys.stderr)
	assert log.isEnabledFor(INFO)
