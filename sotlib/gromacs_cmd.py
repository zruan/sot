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

'''Utility functions to call gromacs/plumed command
'''

import logging
import subprocess
from os import path
from glob import glob

def pdb2gmx(pdb_image, force_field='amber99sb-ildn', water_model='tip3p', \
            ignh=False, vsites='hydrogen', gromacs='', \
            stdout=subprocess.STDOUT):
    # Wrapper for pdb2gmx command
    assert path.isfile(pdb_image), "Input file (%s) not found!" % pdb_image
    pdb2gmx_cmd = [path.join(gromacs, 'gmx'), 'pdb2gmx',
                   '-f', pdb_image,
                   '-vsite', vsites,
                   '-ff', force_field,
                   '-water', water_model,
                   '-o', pdb_image+'_pdb2gmx.pdb',
                   '-p', pdb_image+'_topol.top',
                   '-i', pdb_image+'_posre.itp']
    # be careful of using the -ignh. In some cases, the histidine residue
    # will be protonated differently along the string images, causing
    # indexing issues.
    # TODO: disable the option?
    if ignh: pdb2gmx_cmd.append('-ignh')
    logging.info('running pdb2gmx: ' + ' '.join(pdb2gmx_cmd))
    subprocess.check_call(pdb2gmx_cmd, stdout=stdout, \
                          stderr=subprocess.STDOUT)

def editconf(pdb_image, distance='1.0', box_type='dodecahedron', center=True, \
             gromacs='', stdout=subprocess.STDOUT):
    # Wrapper for editconf command
    assert path.isfile(pdb_image+'_pdb2gmx.pdb'), \
            "Input file (%s) not found!" % pdb_image+'_pdb2gmx.pdb'
    editconf_cmd = [path.join(gromacs, 'gmx'), 'editconf',
                   '-f', pdb_image+'_pdb2gmx.pdb',
                   '-o', pdb_image+'_box.pdb',
                   '-d', distance,
                   '-bt', box_type]
    if center: editconf_cmd.append('-c')
    logging.info('running editconf: ' + ' '.join(editconf_cmd))
    subprocess.check_call(editconf_cmd, stdout=stdout, \
                              stderr=subprocess.STDOUT)

def solvate(pdb_image, library='spc216.gro', gromacs='', \
            stdout=subprocess.STDOUT):
    # Wrapper for solvate command (Gromacs 5.0+)
    assert path.isfile(pdb_image+'_box.pdb'), \
            "Input file (%s) not found!" % pdb_image+'_box.pdb'
    genbox_cmd = [path.join(gromacs, 'gmx'), 'solvate',
                  '-cp', pdb_image+'_box.pdb', 
                  '-cs', 'spc216.gro',
                  '-o', pdb_image+'_solv.pdb', 
                  '-p', pdb_image+'_topol.top']
    logging.info('solvation: ' + ' '.join(genbox_cmd))
    subprocess.check_call(genbox_cmd, stdout=stdout, \
                              stderr=subprocess.STDOUT)
def grompp(pdb_image, mdp, coordinate, topology, outsuffix, gromacs='', \
           stdout=subprocess.STDOUT):
    # Wrapper for grompp command
    assert path.isfile(mdp), "Input mdp file (%s) not found!" % mdp
    assert path.isfile(coordinate), \
            "Input coordinate file (%s) not found!" % coordinate
    grompp_cmd = [path.join(gromacs, 'gmx'), 'grompp',
                  '-f', mdp, 
                  '-c', coordinate,
                  '-p', topology,
                  '-o', pdb_image+outsuffix+'.tpr',
                  '-po', pdb_image+outsuffix+'.mdp']
    logging.info('running grompp: ' + ' '.join(grompp_cmd))
    subprocess.check_call(grompp_cmd, stdout=stdout, \
                              stderr=subprocess.STDOUT)

def genion(pdb_image, conc='0.1', input='13', gromacs='', \
        stdout=subprocess.STDOUT):
    # Wrapper for genion command
    assert path.isfile(pdb_image+'_ions.tpr'), \
            "Input file (%s) not found!" % (pdb_image+'_ions.tpr')
    genion_cmd = [path.join(gromacs, 'gmx'), 'genion',
                  '-s', pdb_image+'_ions.tpr',
                  '-o', pdb_image+'_ions.pdb', 
                  '-p', pdb_image+'_topol.top',
                  '-conc', conc,
                  '-neutral']
    logging.info('adding ions: ' + ' '.join(genion_cmd))
    p = subprocess.Popen(genion_cmd, stdin=subprocess.PIPE, \
                             stdout=stdout, stderr=subprocess.STDOUT)
    p.communicate(input=input)

def mdrun(pdb_image, suffix, gromacs='', plumed=None, gpuid=None, nt=None, pin=None, \
        stdout=subprocess.STDOUT):
    # Wrapper for mdrun command
    mdrun_cmd = [path.join(gromacs, 'gmx'), 'mdrun',
                 '-s', pdb_image+suffix+'.tpr',
                 '-o', pdb_image+suffix+'.trr',
                 '-x', pdb_image+suffix+'.xtc',
                 '-c', pdb_image+suffix+'.pdb',
                 '-g', pdb_image+suffix+'.log',
                 '-cpo', pdb_image+suffix+'.cpt',
                 '-e', pdb_image+suffix+'.edr']
    if plumed:
        mdrun_cmd.extend(['-plumed', plumed])
    if gpuid:
        mdrun_cmd.extend(['-gpu_id', gpuid])
    if nt:
        mdrun_cmd.extend(['-nt', str(nt)])
    if pin:
        mdrun_cmd.extend(['-pin', 'on'])
    logging.info('mdrun command: ' + ' '.join(mdrun_cmd))
    subprocess.check_call(mdrun_cmd, stdout=stdout, \
                              stderr=subprocess.STDOUT)

def make_ndx(base, file_identifier, input='q\n', gromacs='', \
        stdout=subprocess.STDOUT):
    # Wrapper for make_ndx command
    tpr_files = glob(path.join(base, '*'+file_identifier+'*.tpr'))
    assert len(tpr_files) != 0, "No tpr files found (%s)!" % \
            path.join(base, '*'+file_identifier+'*.tpr')
    # we just need one tpr file to create the index
    make_ndx_cmd = [path.join(gromacs, 'gmx'), 'make_ndx',
                    '-f', tpr_files[0],
                    '-o', path.join(base, 'index.ndx')]
    logging.info('writing index.ndx: ' + ' '.join(make_ndx_cmd))
    #p = subprocess.Popen(make_ndx_cmd, stdin=subprocess.PIPE)
    p = subprocess.Popen(make_ndx_cmd, stdin=subprocess.PIPE, \
                             stdout=stdout, stderr=subprocess.STDOUT)
    p.communicate(input=input)

def trjcat(base, file_identifier, input='1', trj_suffix='.xtc', \
        index='index.ndx', out_file='swarm.xtc', gromacs='', \
        stdout=subprocess.STDOUT):
    # Wrapper for trjcat command
    trj_files = glob(path.join(base, '*'+file_identifier+'*'+trj_suffix))
    assert len(trj_files) != 0, "No trajectory files found (%s)!" % \
            path.join(base, '*'+file_identifier+'*'+trj_suffix)
    assert path.isfile(path.join(base, index)), \
            "Input file (%s) not found!" % path.join(base, index)
    # Now concatenate the trajectories
    trjcat_cmd = [path.join(gromacs, 'gmx'), 'trjcat',
                    '-cat',
                    '-n', path.join(base, index),
                    '-o', path.join(base, out_file),
                    '-f']
    trjcat_cmd.extend(trj_files)
    logging.info('Concatenating trajectories: ' + ' '.join(trjcat_cmd))
    p = subprocess.Popen(trjcat_cmd, stdin=subprocess.PIPE, \
                             stdout=stdout, stderr=subprocess.STDOUT)
    p.communicate(input=input)


def plumed(base, trajectory, plumed, directory='', format='xtc', stdout=subprocess.STDOUT):
    # Wrapper for plumed command
    plumed_cmd = [path.join(directory, 'plumed'), 'driver',
                 '--plumed', plumed]
    if format == 'xtc':
        plumed_cmd.extend(['--mf_xtc', path.join(base, trajectory)])
    elif format == 'pdb':
        plumed_cmd.extend(['--mf_pdb', path.join(base, trajectory)])
    logging.info('plumed command: ' + ' '.join(plumed_cmd))
    subprocess.check_call(plumed_cmd, stdout=stdout, \
                              stderr=subprocess.STDOUT)

def cat(file_list, output):
    cmd = 'cat ' + ' '.join(file_list) + ' > ' + output 
    logging.info('concatenating pdb files: ' + cmd)
    subprocess.check_call(cmd, shell=True)


def trjconv(base, trj_file, out_file, top_file=None, inputc=None, idx=None, \
        args=None, gromacs='', stdout=subprocess.STDOUT):
    # Wrapper for trjconv command
    trj_files = glob(path.join(base, trj_file))
    assert path.isfile(path.join(base, trj_file)), \
            "Input trajecotry file (%s) not found!" % path.join(base, trj_file)
    # Now concatenate the trajectories
    trjconv_cmd = [path.join(gromacs, 'gmx'), 'trjconv',
                    '-f', path.join(base, trj_file),
                    '-o', path.join(base, out_file)]
    if top_file:
        trjconv_cmd.extend(['-s', path.join(base, top_file)])
    if idx:
        assert path.isfile(path.join(base, idx)), \
                "Input index file (%s) not found!" % path.join(base, idx)
        trjconv_cmd.extend(['-n', path.join(base, idx)])
    if args:
        trjconv_cmd.extend(args)
    logging.info('Converting trajectories: ' + ' '.join(trjconv_cmd))
    p = subprocess.Popen(trjconv_cmd, stdin=subprocess.PIPE, \
                             stdout=stdout, stderr=subprocess.STDOUT)
    p.communicate(input=inputc)


