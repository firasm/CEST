version = '1.0'
release = version + 'alpha'
from sarpy.helpers import git_repo_state
repo_state = git_repo_state()

import logging
logger=logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())