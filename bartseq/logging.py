from logging import (
	Logger, basicConfig, lastResort,
	DEBUG, INFO, WARNING, ERROR, CRITICAL,
)


log = Logger(__name__)


def init_logging():
	basicConfig(level=INFO)
	lastResort.setLevel(INFO)
