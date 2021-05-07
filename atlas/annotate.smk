import os
import re
import sys
import tempfile
#import pandas as pd
#import numpy as np

from snakemake.utils import logger, min_version

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"scripts"))
import utils

from conf import update_config
config= update_config(config)

TMPDIR = config.get('tmpdir',tempfile.gettempdir() )

#CONDAENV = "envs" # overwrite definition in download.smk

include: "rules/dram.smk"
