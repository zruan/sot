#!/usr/bin/env python

###############################################################################
# Swarms of trajectory (SOT) -- https://github.com/zruan/sot                  #
#                                                                             #
# Copyright (c) 2017 Zheng Ruan and the University of Georgia. All rights     #
# reserved.                                                                   #
#                                                                             #
# This program is free software; you can redistribute it and/or               #
# modify it under the terms of the GNU General Public License                 #
# as published by the Free Software Foundation; either version 2              #
# of the License, or (at your option) any later version.                      #
#                                                                             #
###############################################################################

"""SOT -- A python implementation of the string method with swarms of
trajectories using GROMACS as backend molecular dynamics engine
"""

from os.path import dirname, join

DIR = (dirname(__file__) or '.')

setup_args = dict(
    name='sot',
    version='0.1',
    description=__doc__,
    author='Zheng Ruan',
    author_email='zruan1991@gmail.com',
    url='https://github.com/zruan/sot',
    classifiers=[
        'Development Status :: Pre - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GPL License 2.0',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Computational Chemistry',
    ],
    packages=['sotlib'],
)

try:
    from setuptools import setup
    setup_args.update(
        install_requires=[
            'mdtraj >= 1.9.0',
        ])
except:
    from distutils.core import setup

setup(**setup_args)

