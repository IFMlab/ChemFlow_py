'''
This file is required to run the gromacs module
To get this script follow the installation instructions and requirements for gromacs.

You can request this file from: http://mackerell.umaryland.edu/CHARMM_ff_params.html

Here is the header of the file, if you want to check the compatibility:
'''

#!/usr/bin/env python3

# USAGE: python cgenff_charmm2gmx.py DRUG drug.mol2 drug.str charmm36.ff
# Tested with Python 3.5.2. Requires numpy and networkx
# The networkx version MUST be in the 1.x series. Tested version: 1.11

# Copyright (C) 2014 E. Prabhu Raman prabhu@outerbanks.umaryland.edu
#
# Modified 11/6/2018 by Justin Lemkul to add lone pair support
# needed for CGenFF >= 4.0 halogens
#
# Modified 01/10/2019 by Conrard Tetsassi to work with Networkx 2.3
# Included notes on bonds, angles, and dihedrals
#
# For help/questions/bug reports, please contact Justin Lemkul jalemkul@vt.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU Affero General Public License for more details.
# <http://www.gnu.org/licenses/>

# EXAMPLE: You have a drug-like molecule in drug.mol2 file
# ParamChem returns a CHARMM stream file drug.str with topology and parameters
# INPUT
# The program needs four inputs:
#	(1) The first argument (resname) is found in the RESI entry in the CHARMM stream file; for example
#		RESI DRUG		  0.000  ! here DRUG is the resname
#	(2) drug.mol2 is the .mol2 which you supplied to ParamChem
#	(3) drug.str contains the topology entry and CGenFF parameters that ParamChem generates for you
#	(4) charmm36.ff should contain the CHARMM force-field in GROMACS format
#		Download it from: http://mackerell.umaryland.edu/CHARMM_ff_params.html

# OUTPUT
# The program will generate 4 output files ("DRUG" is converted to lowercase and the files are named accordingly):
#	(1) drug.itp - contains GROMACS itp
#	(2) drug.prm - contains parameters obtained from drug.str which are converted to GROMACS format and units
#	(3) drug.top - A Gromacs topology file which incorporates (1) and (2)
#	(4) drug_ini.pdb - Coordinates of the molecule obtained from drug.mol2

# The program has been tested only on CHARMM stream files containing topology and parameters of a single molecule.