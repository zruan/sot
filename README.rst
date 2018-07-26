============================
Swarms of trajectories (SOT)
============================

A python implementation of the string method with swarms-of-trajectories (SOT)
with GROMACS as the backend molecular dynamics engine.

Freely distributed under the GNU General Public License v2 (see LICENSE).

The code is not estensively tested. Make sure you understand the code well before
carrying out productive simulations.

Installation
------------

The installation consists of a Python library, ``sotlib``.

Download the .zip file and unpack it, or clone this git repository, to get the
source code.

To use all the features of sot, you'll need the following third-party
programs installed:

- Python_ 2.7
- GROMACS_ 5.0+
- PLUMED_ 2.2+
- MDtraj_ 1.9.0+

.. _Python: http://www.python.org/download/
.. _GROMACS: http://www.gromacs.org/Downloads
.. _PLUMED: http://www.plumed.org/get-it
.. _MDtraj: http://mdtraj.org/1.9.0/

Note that PLUMED must be incorporated into GROMACS through patching. Users are
encouraged to read the online `PLUMED documentation
<https://plumed.github.io/doc-v2.3/user-doc/html/_installation.html#InstallingPlumed>`_.

Once the above programs are installed, run the setup file in SOT root folder:

    python setup.py build
    python setup.py install

Basic usage
-----------

Global options:

  ``-h``, ``--help``
      Show a help message and basic usage.

Sub-commands:

    `minimize`_
        Minimize the input string images
    `thermo`_
        Thermalization of the system under specific thermostat
    `equil`_
        Equilibrate the system under specific pressure
    `swarm`_
        Generate short swarm of trajectories
    `update-string`_
        Generate a new string based on the average CVs from
        swarm runs.

Commands
--------

minimize
````````

Minimize a list of input string images.

Assume we have a directory ``iter0`` containing a set of string images
(`path_0.pdb`, `path_1.pdb`, ...). The ``minimize`` sub-command calls
`pdb2gmx`, `editconf`, `solvate`, `grompp`, `genion`, `grompp`, and `mdrun`.
from ``GROMACS`` to prepare an MD system and minimize the system energy.

A few important options shall be specified in the ``minimize`` step:

    **--mdp-minim**
        The mdp configuration file for gromacs to run minimization
    **--mdp-minim2**
        The mdp configuration file for the second round of gromacs minimization

In most case, I find a single steepest descent minimization is enough, since we
are doing swarms of short simulation. However, if the user prefers, a conjugate
gradient minimization can be carried out.

The ``minimize`` sub-command also allows the user to choose force field and
water models for the simulation. This is achieved by specifying the following
options:

    **--ff**
        The force field provided by GROMACS. The default force field is
        `amber99sb-ildn`
    **--water**
        The water model for the simulation. Default one is ``tip3p``

If a more advanced control over the MD system construction is needed. The user
shall hack into the code. Specifically the `func_minimize` defined in
``sot/sotlib/swarmlib.py`` is the main routine that carry out the minimization.
In addition, the ``sot/sotlib/gromacs_cmd.py`` calls the GROMACS program by
passing the user specified options. Advanced users may also want to modify the
way in which GROMACS commands get called.

Finally, the **--image-prefix** specifies the string image prefix. In this
example, the prefix shall be `path_`. The program will automatically detects
files that match `path_[0-9]+.pdb` and extract the number associated with them.

An example command can be::
    
    python /path/to/swarm.py minimize --mdp-minim mdp/minim.mdp --image-prefix path_ --ff amber99sb-ildn --water tip3p --clear iter0/

The command will run energy minimization and generate `path_[0-9].pdb_em.pdb`
file for each string images under ``iter0/``. These are the minimized pdb
structures. The **--clear** option will remove the `.trr` file from the base
directory as they are space consuming.

thermo
``````

Bring temperature to the system by coupling to a specific thermostat.

This is typically the following step for the energy minimized string images.
The only important option at this step is the mdp file for thermolization
(specified by **--mdp-thermo**). Note that pressure coupling should be turned
off in this step.

An example command can be::
    
    python /path/to/swarm.py thermo --image-prefix path_ --mdp-thermo mdp/nvt.mdp --clear iter0/

The command will run thermolization and generate `path_[0-9]_nvt.pdb` file
for each string images under ``iter0/``. These are the thermolized pdb
structures.

equil
`````

Run restrained sampling under specific barostat.

This is typically the following step for the thermolized string images.
The only important option at this step is the mdp file for equilibration
(specified by **--mdp-equil**).

An example command can be::
    
    python /path/to/swarm.py thermo --image-prefix path_ --mdp-equil mdp/npt.mdp --clear iter0/

The command will run equilibration and generate `path_[0-9]_npt.pdb` file
for each string images under ``iter0/``. These are the equilibrated pdb
structures.

swarm
`````

The `swarm` sub-command launches a number of short unbiased trajectories, each
starting at the equilibrated string images. The CVs shall be defined in a
`plumed.dat` file and can be calculated on-the-fly through **--plumed** option.
Several important options are explained below:

    **--num**
        The number of short swarm trajectories to run
    **--plumed**
        plumed configuration file that defines the CVs. Note that `RESTART`
        keyword should not be enabled in the file.
    **--plumed-output**
        The output file defined in the plumed configuration file. The file name
        should match exactly.
    **--mdp-swarm**
        The GROMACS mdp file to run unbiased swarm MDs. Since a large number of
        swarm trajectories will be running, it is recommended to output
        trajectory in `.xtc` format. To achieve it, the `nstxout`, `nstvout`,
        `nstenergy`, `nstlog` should all be 0. The frequency of trajectory
        output in `.xtc` file is controled solely by `nstxtcout`. In addition,
        the velocity generation shall be turned on (`gen-val=yes`) with random
        seed (`gen-seed`) set to -1. Refer `GROMACS documentation for further
        details<http://manual.gromacs.org/online/mdp_opt.html>`_.

An example command can be::
    
    python /path/to/swarm.py swarm --image-prefix path_ --mdp-swarm mdp/swarm.mdp --num 20 --plumed plumed.dat --plumed-output iter0/COLVAR iter0/

The command will run 20 short unbiased MD simulations each initiated with
different starting velocity. The CVs will also be calculated and averaged for
all the snapshots to update the string images. A `.npy` file will be written
in the base directory, storing the updated CVs.

update-string
`````````````

Find the closest snapshots to the update strings.

Specifically, the sub-command will read the new string images defined from the
SOT simulation and identify the shapshots that are closest to each string images
from the SOT simulation. The identified snapshots will be saved into a new directory
for the next iteration.

