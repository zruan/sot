#!/usr/bin/python2.7

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

"""The package implement the string method with swarm of trajectories by
Benoit Roux [1]_

The script accepts an initial string over the collective variable (CV) space.
And iteratively performs:

    - Energy minimization of the images along the initial string
    - Restrained thermalization of each image (NVT ensemble)
    - Equilibration of each image in NPT ensemble
    - Run a number of short swarm of trajecotires for each image using
         different random seeds
    - Collect statistics of CVs from swarm simulations and Update the string
         for subsequent iteraction


.. [1] J. Phys. Chem. B, 2008, 112 (11), pp 3432-3440

"""

import logging
import argparse

from sotlib.swarmlib import func_minimize, func_thermo, func_equilibration, func_swarm, func_update_string

if __name__ == "__main__":
    AP = argparse.ArgumentParser(
            description=__doc__,
            epilog="Contact Zheng Ruan <rzzmh@uga.edu> for help.")

    # Global options
    AP.add_argument('-q', '--quiet',
            action='store_true',
            help="Don't print status messages, only warnings and errors.")
    AP_subparsers = AP.add_subparsers(
            help='sub-commands (use with -h for more info)')

    # Sub-command: func_minimize
    P_minimize = AP_subparsers.add_parser('minimize',
                help='Minimize the input string images')
    P_minimize.add_argument('basedir',
                help='Top of the directory tree where the string images are')
    P_minimize.add_argument('--image-prefix', dest='image_prefix',
                required=True,
                help="Prefix of string images. Input images has to be in pdb format")
    P_minimize.add_argument('--mdp-minim', dest='mdp_minim',
                required=True, help='Gromacs mdp file for energy minimization')
    P_minimize.add_argument('--mdp-minim2', dest='mdp_minim2', default=None,
                help='Gromacs mdp file for second round of minimization energy')
    P_minimize.add_argument('--ff', dest='force_field',
                default='amber99sb-ildn', help='Force field for MD.')
    P_minimize.add_argument('--water', dest='water', default='tip3p',
                help="Water models to use")
    P_minimize.add_argument('--gromacs-bin', dest='gromacs_bin', default='',
                help="Directory containing binary files for gromacs")
    P_minimize.add_argument('--clear', dest='clear', action='store_true',
                default=False, help="Weather to clear up the output")
    P_minimize.add_argument('--nt', dest='nt', default=None,
                help="Number of cores to use for mdrun")
    P_minimize.add_argument('--pin', dest='pin', action='store_true',
                help="Weather to pin threads")
    P_minimize.add_argument('--gpuid', dest='gpuid', default=None,
                help="The GPU ID for mdrun")
    P_minimize.set_defaults(func=func_minimize)

    # Sub-command: thermalization
    P_thermo = AP_subparsers.add_parser('thermo',
                help='Thermalization of the system under specific thermostat')
    P_thermo.add_argument('basedir',
                help='Top of the directory tree where the string images are')
    P_thermo.add_argument('--image-prefix', dest='image_prefix',
                required=True, help="Prefix of energy minimized string images.\
                                     Two files are expected: \
                                     image_prefix[0-9]+.*[gro|pdb] and \
                                     image_prefix[0-9]+.*top")
    P_thermo.add_argument('--mdp-thermo', dest='mdp_thermo', required=True,
                help='Gromacs mdp file for thermalization')
    P_thermo.add_argument('--struc-suffix', dest='struc_suffix',
                default='.pdb_em', help='The suffix of input structures')
    P_thermo.add_argument('--continue', dest='cont', action='store_true',
                help='Run in continue mode')
    P_thermo.add_argument('--gromacs-bin', dest='gromacs_bin', default='',
                help="Directory containing binary files for gromacs")
    P_thermo.add_argument('--clear', dest='clear', action='store_true',
                default=False, help="Weather to clear up the output")
    P_thermo.add_argument('--nt', dest='nt', default=None,
                help="Number of cores to use for mdrun")
    P_thermo.add_argument('--pin', dest='pin', action='store_true',
                help="Weather to pin threads")
    P_thermo.add_argument('--gpuid', dest='gpuid', default=None,
                help="The GPU ID for mdrun")
    P_thermo.set_defaults(func=func_thermo)

    # Sub-command: equilibration
    P_equil = AP_subparsers.add_parser('equil',
                help='Equilibrate the system under specific pressure')
    P_equil.add_argument('basedir',
                help='Top of the directory tree where the string images are')
    P_equil.add_argument('--image-prefix', dest='image_prefix',
                required=True, help="Prefix of thermalized string images.\
                                     Two files are expected: \
                                     image_prefix[0-9]+.*[gro|pdb] and \
                                     image_prefix[0-9]+.*top")
    P_equil.add_argument('--mdp-equil', dest='mdp_equil', required=True,
                help='Gromacs mdp file for equilibration')
    P_equil.add_argument('--struc-suffix', dest='struc_suffix',
                default='_nvt', help='The suffix of input structures')
    P_equil.add_argument('--continue', dest='cont', action='store_true',
                help='Run in continue mode')
    P_equil.add_argument('--gromacs-bin', dest='gromacs_bin', default='',
                help="Directory containing binary files for gromacs")
    P_equil.add_argument('--clear', dest='clear', action='store_true',
                default=False, help="Weather to clear up the output")
    P_equil.add_argument('--nt', dest='nt', default=None,
                help="Number of cores to use for mdrun")
    P_equil.add_argument('--pin', dest='pin', action='store_true',
                help="Weather to pin threads")
    P_equil.add_argument('--gpuid', dest='gpuid', default=None,
                help="The GPU ID for mdrun")
    P_equil.set_defaults(func=func_equilibration)

    # Sub-command: swarm
    P_swarm = AP_subparsers.add_parser('swarm',
                help='Generate short swarm of trajectories')
    P_swarm.add_argument('basedir',
                help='Top of the directory tree where the string images are')
    P_swarm.add_argument('--image-prefix', dest='image_prefix',
                required=True, help="Prefix of equilibrated string images.\
                                     Two files are expected: \
                                     image_prefix[0-9]+.*[gro|pdb] and \
                                     image_prefix[0-9]+.*top")
    P_swarm.add_argument('--num', dest='num', required=True, type=int,
                help='Number of swarm simulations to run')
    P_swarm.add_argument('--plumed', dest='plumed', required=True,
                help='plumed file that defines the CVs (2 CVs expected)')
    P_swarm.add_argument('--plumed-output', dest='plumed_output',
                default='COLVAR',
                help='plumed file that defines the CVs (2 CVs expected)')
    P_swarm.add_argument('--mdp-swarm', dest='mdp_swarm', required=True,
                help='Gromacs mdp file for swarm simulation')
    P_swarm.add_argument('--struc-suffix', dest='struc_suffix',
                default='_npt', help='The suffix of input structures')
    P_swarm.add_argument('--continue', dest='cont', action='store_true',
                help='Run in continue mode')
    P_swarm.add_argument('--gromacs-bin', dest='gromacs_bin',
                default='',
                help="Directory containing binary files for gromacs")
    P_swarm.add_argument('--output-suffix', dest='output_suffix', default='string',
                help="File suffix to store the new string CVs (.npy file)")
    P_swarm.add_argument('--nt', dest='nt', default=None,
                help="Number of cores to use for mdrun")
    P_swarm.add_argument('--pin', dest='pin', action='store_true',
                help="Weather to pin threads")
    P_swarm.add_argument('--gpuid', dest='gpuid', default=None,
                help="The GPU ID for mdrun")
    # Do to allow trajectory clearance in this step (do this in the next step)
    #P_swarm.add_argument('--clear', dest='clear', action='store_true',
    #            default=False, help="Weather to clear up the .trr file")
    P_swarm.set_defaults(func=func_swarm)

    # Sub-command: update-string
    P_update = AP_subparsers.add_parser('update-string', 
                help="""Generate a new string based on the average CVs from swarm
                runs.""")
    P_update.add_argument('basedir',
                help='Top of the directory tree where the string images are')
    P_update.add_argument('--newdir', dest='newdir', required=True,
                help='Directory to store the updated string of images')
    P_update.add_argument('--new-image-prefix', dest='new_image_prefix',
                required=True,
                help="Prefix of string images. Input images has to be in pdb format")
    P_update.add_argument('--cvfile', dest='cvfile', required=True,
                help='.npy file that stores the CV values from swarm of trajectories')
    P_update.add_argument('--plumed', dest='plumed', required=True,
                help='Plumed file that defines the CVs (2 CVs expected)')
    P_update.add_argument('--plumed-output', dest='plumed_output',
                default='COLVAR',
                help='The argument must match the one defined in plumed configuration file')
    P_update.add_argument('--continue', dest='cont', action='store_true',
                help='Run in continue mode')
    P_update.add_argument('--gromacs-bin', dest='gromacs_bin', default='',
                help="Directory containing binary files for gromacs")
    P_update.add_argument('--swarm', dest='swarm', default='swarm',
                help='File identifier for the swarm trajectories')
    P_update.add_argument('--clear', dest='clear', action='store_true',
                default=False, help="Weather to clear up the trajectory file")
    P_update.set_defaults(func=func_update_string)

    args = AP.parse_args()
    # Handle global options here
    if args.quiet:
        logging.basicConfig(level=logging.WARNING,
                format="%(module)s: %(message)s")
    else:
        logging.basicConfig(level=logging.INFO,
                format="%(module)s [@%(lineno)s]: %(message)s")
    # Pass the rest along to the sub-command implementation
    args.func(args)


