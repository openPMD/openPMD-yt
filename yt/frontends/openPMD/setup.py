#!/usr/bin/env python
"""

Copyright (c) 2015, Daniel Grassinger (HZDR)

"""

import setuptools
import os
import sys
import os.path


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('openPMD', parent_package, top_path)
    config.make_config_py()  # installs __config__.py
    #config.make_svn_version_py()
    return config
