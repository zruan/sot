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

'''Sub-commands for SOT analysis
'''

from __future__ import print_function
import sys
import logging
import re
import fileinput
import numpy as np
import functools
import shutil
import matplotlib.pyplot as plt
import copy

from os import path, devnull, remove, mkdir
from glob import glob
from itertools import combinations

from gromacs_cmd import *

from reparametrize import rep_pts

def string_dist(string):
    norm = np.max(string, axis=0) - np.min(string, axis=0)
    length = string.shape[0]
    dist = 0 
    for i in range(length-1):
        diff = string[i,:] - string[i+1,:]
        dist += ((diff[0]/norm[0])**2 + (diff[1]/norm[1])**2)**(0.5)
    return dist


def swap(string, num1, num2):
    tmp = copy.copy(string[num1,:])
    string[num1,:]=string[num2,:]
    string[num2,:]=tmp
    return string


def reorder_string(string):
    length = string.shape[0]
    num_pair = combinations(range(0, length), 2)
    for num in num_pair:
        tmp_string = copy.copy(string)
        swap_string = swap(tmp_string, num[0], num[1])
        if string_dist(string) > string_dist(swap_string):
            return swap_string
    return string


def repeat_reorder_string(string, times=10):
    updated_string = [string]
    while times != 0:
        new_string = reorder_string(updated_string[-1])
        updated_string.append(new_string)
        times -= 1
    return updated_string[-1]


def remove_files(basedir, name, suffix='.trr'):
    fname = []
    for f in glob(path.join(basedir, '*'+name+'*'+suffix)):
        fname.append(f)
        remove(f)
    return fname


def check_basedir(args):
    assert path.isdir(args.basedir.rstrip('/'))
    return args.basedir.rstrip('/')


def check_im_consecutivity(string_images):
    numbers = string_images.keys()
    assert set(numbers) == set(range(min(numbers), max(numbers)+1)), \
            'String image numbers are not consecutive. %s detected.' % str(numbers)
    return numbers


def check_convergence(logfile):
    with open(logfile, 'r') as logf:
        for line in logf.readlines():
            if 'converged to Fmax <' in line:
                return True
    return False


def func_minimize(args):
    '''Function to perform energy minimization along a string.

    The function execude a sequential of gromacs command to prepare the tpr
    file and use mdrun to minimize the system energy.
    '''
    logging.info('### Start func_minimize ###')
    base = check_basedir(args)

    # Identify string images (in pdb format) from the base folder
    # String images have to be sequentially numbered in the basedir
    # e.g.: image_0.pdb, image_1.pdb, image_2.pdb, ...
    # with 'image_' as image_prefix
    image_prefix = args.image_prefix
    string_re = re.compile('^'+image_prefix+'(?P<num>[0-9]+).pdb$')
    string_images = {}
    for i in glob(path.join(base, image_prefix)+'*pdb'):
        match = string_re.match(path.basename(i))
        if match:
            string_images[int(match.groupdict()['num'])] = i

    # check consecutivity of image number 
    numbers = check_im_consecutivity(string_images)
    logging.info('%d string images processed.' % len(numbers))

    # Process string image files using gromacs
    gromacs_bin = args.gromacs_bin
    force_field = args.force_field
    water_model = args.water
    mdp_minim = args.mdp_minim
    FNULL = open(devnull, 'w')
    for num, image in string_images.items():
        logging.info('Processing image %d (%s)' % (num, image))
        # pdb2gmx
        pdb2gmx(image, force_field=force_field, water_model=water_model, \
                ignh=True, vsites='hydrogen', gromacs=gromacs_bin, stdout=FNULL)
        # editconf
        editconf(image, distance='1.0', box_type='dodecahedron', center=True, \
                 gromacs=gromacs_bin, stdout=FNULL)
        # solvate (gromacs 5.0+)
        solvate(image, library='spc216.gro', gromacs='', stdout=FNULL)
        # grompp 
        grompp(image, mdp=mdp_minim, coordinate=image+'_solv.pdb', \
               topology=image+'_topol.top', outsuffix='_ions', \
               gromacs=gromacs_bin, stdout=FNULL)
        # genion
        genion(image, conc='0.1', gromacs=gromacs_bin, stdout=FNULL)
        # grompp 
        grompp(image, mdp=mdp_minim, coordinate=image+'_ions.pdb', \
               topology=image+'_topol.top', outsuffix='_em', \
               gromacs=gromacs_bin, stdout=FNULL)
        # mdrun
        mdrun(image, suffix='_em', nt=args.nt, pin=args.pin, gpuid=args.gpuid, \
               gromacs=gromacs_bin, stdout=FNULL)

        # .top file will include the .itp file from basedir/
        # This will cause problem for the subsequent grompp command
        # Here is the dirty fix
        if path.isfile(image+'_topol.top'):
            for line in fileinput.input([image+'_topol.top'], inplace=True, backup='.bak'):
                # An example record: include "iter0/path_0.pdb_posre.itp"
                if base+'/' in line:
                    print(line.replace(base+'/', '').rstrip('\n'))
                else:
                    print(line.rstrip('\n'))

        # check to see convergence
        if check_convergence(image+'_em.log'):
            logging.info('Energy minimization converged (%s)' % (image+'_em.log'))
            if args.clear:
                fn = remove_files(base, 'em')
                logging.info('Files removed (%s)' % ','.join(fn))
        else:
            if args.mdp_minim2:
                # grompp 
                grompp(image, mdp=args.mdp_minim2, coordinate=image+'_em.pdb', \
                       topology=image+'_topol.top', outsuffix='_em2', \
                       gromacs=gromacs_bin, stdout=FNULL)
                # mdrun
                mdrun(image, suffix='_em2', nt=args.nt, pin=args.pin, gpuid=args.gpuid, \
                       gromacs=gromacs_bin, stdout=FNULL)
                # dirty solution: we just copy the output of _em2.pdb to _em.pdb
                shutil.copyfile(image+'_em2.pdb', image+'_em.pdb')
                if check_convergence(image+'_em2.log'):
                    logging.info('Energy minimization converged (%s)' % (image+'_em2.log'))
                    if args.clear:
                        fn = remove_files(base, 'em')
                        logging.info('Files removed (%s)' % ','.join(fn))
    return


def func_thermo(args):
    '''Function to perform thermalization along a string.

    The function processes energy minimized structures and use mdrun to bring
    temperatue to the system.
    '''
    logging.info('### Start func_thermo ###')
    base = check_basedir(args)

    # Identify minimized topologies and structures
    # topology naming: image_prefix[0-9]+.*[top]
    # structures naming: image_prefix[0-9]+.*[gro|pdb]
    image_prefix = args.image_prefix
    string_topo_re = re.compile('^'+image_prefix+'(?P<num>[0-9]+).*top$')
    string_images = {}
    structure_suffix = args.struc_suffix
    for i in glob(path.join(base, image_prefix)+'*top'):
        match = string_topo_re.match(path.basename(i))
        num = match.groupdict()['num']
        if match:
            struc_file = path.join(base, image_prefix+num+structure_suffix)
            if path.isfile(struc_file+'1.pdb'):
                string_images[int(num)] = (i, struc_file+'1.pdb')
            elif path.isfile(struc_file+'.pdb'):
                string_images[int(num)] = (i, struc_file+'.pdb')
        else:
            raise IOError("%s does't have a corresponding structures (%s expected)" % \
                               (i, struc_file+'.pdb'))

    # check consecutivity of image number 
    numbers = check_im_consecutivity(string_images)
    logging.info('%d minimized string images processed.' % len(numbers))

    # Process string image files using gromacs
    gromacs_bin = args.gromacs_bin
    mdp_thermo = args.mdp_thermo
    FNULL = open(devnull, 'w')
    for num, image in string_images.items():
        logging.info('Processing image %d (%s)' % (num, image))
        # grompp 
        if args.cont and path.isfile(path.join(base, image_prefix+str(num))+'_nvt.pdb'):
            logging.info('Finding trajectory %s, skipping...' % (image_prefix+str(num)+'_nvt.pdb'))
            continue
        else:
            grompp(path.join(base, image_prefix+str(num)), mdp=mdp_thermo, \
                   coordinate=image[1], topology=image[0], outsuffix='_nvt', \
                   gromacs=gromacs_bin, stdout=FNULL)
            # mdrun
            mdrun(path.join(base, image_prefix+str(num)), suffix='_nvt', \
                   nt=args.nt, pin=args.pin, gpuid=args.gpuid, \
                   gromacs=gromacs_bin, stdout=FNULL)
            if args.clear:
                fn = remove_files(base, 'nvt')
                logging.info('Files removed (%s)' % ','.join(fn))
    return


def func_equilibration(args):
    '''Function to perform restrained equilibration along a string.

    The function processes thermalized structures and use mdrun to 
    equilibrate the system under pressure.
    '''
    logging.info('### Start func_equilibration ###')
    base = check_basedir(args)

    # Identify thermalized topologies and structures
    # topology naming: image_prefix[0-9]+.*[top]
    # structures naming: image_prefix[0-9]+.*[gro|pdb]
    image_prefix = args.image_prefix
    string_topo_re = re.compile('^'+image_prefix+'(?P<num>[0-9]+).*top$')
    string_images = {}
    structure_suffix = args.struc_suffix
    for i in glob(path.join(base, image_prefix)+'*top'):
        match = string_topo_re.match(path.basename(i))
        num = match.groupdict()['num']
        if match:
            struc_file = path.join(base, image_prefix+num+structure_suffix)
            if path.isfile(struc_file+'.pdb'):
                string_images[int(num)] = (i, struc_file+'.pdb')
            elif path.isfile(struc_file+'.gro'):
                string_images[int(num)] = (i, struc_file+'.gro')
        else:
            raise IOError("%s does't have a corresponding structures (%s or %s expected)" % \
                               (i, struc_file+'.pdb', struc_file+'.gro'))

    # check consecutivity of image number 
    numbers = check_im_consecutivity(string_images)
    logging.info('%d thermalized string images processed.' % len(numbers))

    # Process string image files using gromacs
    gromacs_bin = args.gromacs_bin
    mdp_equil = args.mdp_equil
    FNULL = open(devnull, 'w')
    for num, image in string_images.items():
        logging.info('Processing image %d (%s)' % (num, image))
        if args.cont and path.isfile(path.join(base, image_prefix+str(num))+'_nvt.pdb'):
            logging.info('Finding trajectory %s, skipping...' % (image_prefix+str(num)+'_nvt.pdb'))
            continue
        else:
            # grompp 
            grompp(path.join(base, image_prefix+str(num)), mdp=mdp_equil, \
                   coordinate=image[1], topology=image[0], outsuffix='_npt', \
                   gromacs=gromacs_bin, stdout=FNULL)
            # mdrun
            mdrun(path.join(base, image_prefix+str(num)), suffix='_npt', \
                  nt=args.nt, pin=args.pin, gpuid=args.gpuid, \
                  gromacs=gromacs_bin, stdout=FNULL)
            if args.clear:
                fn = remove_files(base, 'npt')
                logging.info('Files removed (%s)' % ','.join(fn))
    return


def func_swarm(args):
    '''Function to prepare and run short swarm of trajectories along a string
    of images.
    '''
    logging.info('### Start func_swarm ###')
    base = check_basedir(args)

    # Identify equilibrated topologies and structures
    # topology naming: image_prefix[0-9]+.*[top]
    # structures naming: image_prefix[0-9]+.*[gro|pdb]
    image_prefix = args.image_prefix
    string_topo_re = re.compile('^'+image_prefix+'(?P<num>[0-9]+).*top$')
    string_images = {}
    structure_suffix = args.struc_suffix
    for i in glob(path.join(base, image_prefix)+'*top'):
        match = string_topo_re.match(path.basename(i))
        num = match.groupdict()['num']
        if match:
            struc_file = path.join(base, image_prefix+num+structure_suffix)
            if path.isfile(struc_file+'.pdb'):
                string_images[int(num)] = (i, struc_file+'.pdb')
            elif path.isfile(struc_file+'.gro'):
                string_images[int(num)] = (i, struc_file+'.gro')
        else:
            raise IOError("%s does't have a corresponding structures (%s or %s expected)" % \
                               (i, struc_file+'.pdb', struc_file+'.gro'))

    # check consecutivity of image number 
    numbers = check_im_consecutivity(string_images)
    logging.info('%d equilibrated string images processed.' % len(numbers))

    # Process string image files using gromacs
    num_swarms = args.num
    gromacs_bin = args.gromacs_bin
    mdp_swarm = args.mdp_swarm
    FNULL = open(devnull, 'w')
    cv_center = []
    for num, image in string_images.items():
        cv_mean = []
        for n_swarm in range(num_swarms):
            logging.info('Processing image %d (%s) -- swarm (%d)' % \
                         (num, image, n_swarm))
            # TODO:
            # the args.cont option needs to be rewrite (## not working ##)
            if args.cont and path.isfile(path.join(base, image_prefix+str(num))+'_swarm'+str(n_swarm)+'_protein.pdb'):
                plumed(base, image_prefix+str(num)+'_swarm'+str(n_swarm)+'_protein.pdb',
                       plumed=args.plumed, format='pdb', stdout=FNULL)
            #else:
            else:
                # grompp 
                grompp(path.join(base, image_prefix+str(num)), mdp=mdp_swarm, \
                       coordinate=image[1], topology=image[0], \
                       outsuffix='_swarm'+str(n_swarm), gromacs=gromacs_bin, \
                       stdout=FNULL)
                # mdrun (not using plumed because the hydrogen atoms maybe
                # differently added for each swarms)
                mdrun(path.join(base, image_prefix+str(num)), nt=args.nt,
                       pin=args.pin, gpuid=args.gpuid,
                       suffix='_swarm'+str(n_swarm), gromacs=gromacs_bin,
                       stdout=FNULL)
                # write pdb file for easy processing
                trjconv(base, trj_file=image_prefix+str(num)+'_swarm'+str(n_swarm)+'.xtc',
                       top_file=image_prefix+str(num)+'_swarm'+str(n_swarm)+'.tpr',
                       out_file=image_prefix+str(num)+'_swarm'+str(n_swarm)+'_protein.pdb',
                       inputc='3\n2', # no hydrogen
                       gromacs=gromacs_bin, args=['-fit', 'rot+trans'], stdout=FNULL)
                # remove virtual atoms
                for line in fileinput.input([path.join(base, image_prefix+str(num)+'_swarm'+str(n_swarm)+'_protein.pdb')], \
                                                inplace=True, backup='.bak'):
                    if not (line.startswith('ATOM') and (('MC' in line) or ('MN' in line))):
                        print(line.rstrip('\n'))
                # plumed to calculate CV
                plumed(base, image_prefix+str(num)+'_swarm'+str(n_swarm)+'_protein.pdb',
                       plumed=args.plumed, format='pdb', stdout=FNULL)
            assert path.isfile(args.plumed_output), "CV output file (%s) not found!" % \
                                                     args.plumed_output
            this_cv = np.loadtxt(args.plumed_output)
            cv_mean.append(np.mean(this_cv[:,1:3], axis=0))
            fn = remove_files(base, 'COLVAR', suffix='')
            logging.info('File removed (%s)' % ','.join(fn))
        cv_center.append(np.mean(cv_mean, axis=0))
    np.save('cv_center.npy', np.array(cv_center))
    # reparameterize the string based on the CV center of swarms of trajectories
    reordered_cv_center = repeat_reorder_string(np.array(cv_center))
    result = rep_pts(reordered_cv_center)
    result1 = rep_pts(np.array(result))
    reordered_result = repeat_reorder_string(np.array(result))
    np.save(path.join(base, image_prefix+args.output_suffix+'.npy'), reordered_result)
    return


def func_update_string(args):
    '''Function to extract snapshots based on the given string of images'''

    import mdtraj as md
    from mapping import map_atoms

    logging.info('### Start func_update_string ###')
    base = check_basedir(args)
    assert args.cvfile.endswith('.npy'), \
            'CV file (%s) does not end with .npy!' % args.cvfile
    assert path.isfile(path.join(base, args.cvfile)), \
            'CV file (%s) not found!' % path.join(base, args.cvfile)
    string_cv = np.load(path.join(base, args.cvfile))

    # make ndx file
    gromacs_bin = args.gromacs_bin
    FNULL = open(devnull, 'w')
    ## concatenate trajecory into the single one for later processing
    #if True:
    if not (args.cont and path.isfile(path.join(base, 'together.pdb'))):
        cat(glob(path.join(base, '*swarm*protein*pdb')), path.join(base, 'together.pdb'))
    #    if args.clear:
    #        fn = remove_files(base, args.swarm, suffix='.xtc')
    #        logging.info('Files removed (%s)' % ','.join(fn))
    #    # align trajectory to make molecule whole
        trjconv(base, 'together.pdb', top_file=None, out_file='together.xtc',
                inputc='0', gromacs=gromacs_bin, stdout=FNULL)
    #    # generate 1 pdb file for mdtraj
        top = ''
        with open(path.join(base, 'together.pdb'), 'r') as h:
            for line in h.readlines():
                if 'TER' not in line:
                    top += line
                else:
                    top += line
                    break
        with open(path.join(base, 'together_top.pdb'), 'w') as w:
            w.write(top)
    # run plumed to extract CVs
    plumed(base, 'together.xtc', plumed=args.plumed, stdout=FNULL)

    # identify the closest points in colvar to each string images
    colvar = np.loadtxt(args.plumed_output)
    norm = np.max(colvar[:,1:3], axis=0) - np.min(colvar[:,1:3], axis=0)
    idx = []
    new_string_cv = []
    for image in string_cv:
        dist = ((colvar[:,1]-image[0])/norm[0])**2 + ((colvar[:,2]-image[1])/norm[1])**2
        idx.append((dist).argmin())
        new_string_cv.append(colvar[idx[-1], 1:])
    t = md.load(path.join(base, 'together.xtc'), top=path.join(base, 'together_top.pdb'))
    if not path.isdir(args.newdir): mkdir(args.newdir)
    # remove vsites which will cause problem for next MD iteration
    #vsites_selection = list(t.topology.select('name =~ ".*M[CN].*"'))
    #atoms_to_keep = [a.index for a in t.topology.atoms \
    #                    if a.index not in vsites_selection]
    #t.restrict_atoms(atoms_to_keep)
    for n, i in enumerate(idx):
        logging.info('saving updated string image (%s) to %s.' % \
                         (args.new_image_prefix+str(n)+'.pdb', args.newdir))
        t[i].save(path.join(args.newdir, args.new_image_prefix+str(n)+'.pdb'))
    logging.info('saving CVs of updated string (%s).' % \
                     (path.join(args.newdir, args.new_image_prefix+'new_string.npy')))
    reordered_new_string_cv = repeat_reorder_string(np.array(new_string_cv))
    np.save(path.join(args.newdir, args.new_image_prefix+'new_string.npy'), reordered_new_string_cv)
    # dirty code to copy plumed.dat file
    new_plumed_path = args.plumed.replace(base, args.newdir.rstrip('/'))
    shutil.copyfile(args.plumed, new_plumed_path)
    for line in fileinput.input([new_plumed_path], inplace=True, backup='.bak'):
        # An example record: include "iter0/path_0.pdb_posre.itp"
        if base in line:
            print(line.replace(base, args.newdir.rstrip('/')).rstrip('\n'))
        else:
            print(line.rstrip('\n'))
    #if args.clear:
    #    fn = remove_files(base, 'together_aln.xtc', suffix='')
    #    logging.info('Files removed (%s)' % ','.join(fn))


if __name__ == '__main__':
    pass
