#############################################################################################################################

##       #### ########  ########     ###    ########  #### ########  ######
##        ##  ##     ## ##     ##   ## ##   ##     ##  ##  ##       ##    ##
##        ##  ##     ## ##     ##  ##   ##  ##     ##  ##  ##       ##
##        ##  ########  ########  ##     ## ########   ##  ######    ######
##        ##  ##     ## ##   ##   ######### ##   ##    ##  ##             ##
##        ##  ##     ## ##    ##  ##     ## ##    ##   ##  ##       ##    ##
######## #### ########  ##     ## ##     ## ##     ## #### ########  ######

#############################################################################################################################

import os

import numpy             as np
from pathlib import Path as path
#
### local libraries
#
from ..read.read     import fileio
from ..read.vaspio   import copy_poscar, set_magnetic_moments_from_types, set_magnetic_moments_from_magmom, write_incar, write_potcar
from ..read.vaspio   import get_from_poscar, mag_from_oszcar, mag_from_magcar, mag_from_outcar
from ..thermo.states import permute_spin_staes

#############################################################################################################################

##     ##  #######  ########  ##     ## ##       ########  ######
###   ### ##     ## ##     ## ##     ## ##       ##       ##    ##
#### #### ##     ## ##     ## ##     ## ##       ##       ##
## ### ## ##     ## ##     ## ##     ## ##       ######    ######
##     ## ##     ## ##     ## ##     ## ##       ##             ##
##     ## ##     ## ##     ## ##     ## ##       ##       ##    ##
##     ##  #######  ########   #######  ######## ########  ######

#############################################################################################################################
### Write incar file
"""
potcars = {'path':path/to/paws, 'Cr':'Cr_sv', 'C':'C', 'H':'H'}
types   = {'Cr':{'ls':1.0, 'hs':3.0}} or {'Cr':'Cr3+'}
"""

def vasp(ls:path=None, hs:path=None, ms:list=None, types:list=None, potentials:dict=None) -> None:

    if not isinstance(ms,list):
        if (ms is not None) and (not ms.lower() == "all"): raise RuntimeError("unrecognized option '{ms}'")

    if (ls is None)         and (hs is None): raise RuntimeError("define low- 'ls' and high-spin 'hs' paths")
    if (potentials is None) and (ms is None): raise RuntimeError("'potentials' dictionary not provided")
    if (types is None)      and (ms is None): raise RuntimeError("define ion 'types' as in POSCAR")

    if isinstance(ls,str): ls = path(ls)
    if isinstance(hs,str): hs = path(hs)

    parent_low_spin_state  = path(ls.parent)
    parent_high_spin_state = path(hs.parent)

    if not parent_low_spin_state == parent_high_spin_state: raise RuntimeError(f"parent directory for ls and hs differs")

    if (types is not None) and isinstance(types,list): types = np.asarray(types)

    low_spin_posfile  = ls/"POSCAR"
    low_spin_outfile  = ls/"OUTCAR"
    low_spin_magfile  = ls/"MAGCAR"
    low_spin_oszfile  = ls/"OSZICAR"

    high_spin_posfile = hs/"POSCAR"
    high_spin_outfile = hs/"OUTCAR"
    high_spin_magfile = hs/"MAGCAR"
    high_spin_oszfile = hs/"OSZICAR"
    
    low_spin_poscar   = open(low_spin_posfile).readlines()  if fileio(filename=low_spin_posfile, option="r")  else None
    high_spin_poscar  = open(high_spin_posfile).readlines() if fileio(filename=high_spin_posfile, option="r") else None
 
    low_spin_outcar   = open(low_spin_outfile).readlines()  if fileio(filename=low_spin_outfile, option="s")  else None
    high_spin_outcar  = open(high_spin_outfile).readlines() if fileio(filename=high_spin_outfile, option="s") else None
 
    low_spin_magcar   = open(low_spin_magfile).readlines()  if fileio(filename=low_spin_magfile, option="s")  else None
    high_spin_magcar  = open(high_spin_magfile).readlines() if fileio(filename=high_spin_magfile, option="s") else None
 
    low_spin_oszcar   = open(low_spin_oszfile).readlines()  if fileio(filename=low_spin_oszfile, option="s")  else None
    high_spin_oszcar  = open(high_spin_oszfile).readlines() if fileio(filename=high_spin_oszfile, option="s") else None

    low_spin_atoms, low_spin_elements   = get_from_poscar(data=low_spin_poscar)
    high_spin_atoms, high_spin_elements = get_from_poscar(data=high_spin_poscar)

    if not low_spin_elements   == high_spin_elements:   raise RuntimeError("elements differ")
    if not len(low_spin_atoms) == len(high_spin_atoms): raise RuntimeError("number of atoms differ")

    if types is not None:
        low_spin_magmom = set_magnetic_moments_from_types(spinstate="ls", atoms=low_spin_atoms, elements=low_spin_elements, types=types)

    elif (low_spin_outcar is not None) and (low_spin_oszcar is not None):
        magmom = mag_from_oszcar(data=np.asarray(low_spin_oszcar))
        low_spin_magmom = mag_from_outcar(data=low_spin_outcar, magnetization=magmom)

    elif low_spin_magcar is not None:
        low_spin_magmom = mag_from_magcar(data=low_spin_magcar)

    else:
        raise RuntimeError("insuficient information, define 'types' dictionary, e.g., {'Fe':'Fe3+'}")

    if types is not None:
        high_spin_magmom = set_magnetic_moments_from_types(spinstate="hs", atoms=high_spin_atoms, elements=high_spin_elements, types=types)
    
    elif (high_spin_outcar is not None) and (high_spin_oszcar is not None):
        magmom = mag_from_oszcar(data=np.asarray(high_spin_oszcar))
        high_spin_magmom = mag_from_outcar(data=high_spin_outcar, magnetization=magmom)

    elif high_spin_magcar is not None:
        high_spin_magmom = mag_from_magcar(data=high_spin_magcar)

    else:
        raise RuntimeError("insuficient information, define 'types' dictionary, e.g., {'Fe':'Fe3+'}")

    if ms is None:
        write_incar(filename=ls, magmom=low_spin_magmom)
        write_incar(filename=hs, magmom=high_spin_magmom)

        write_potcar(filename=ls, elements=low_spin_elements,  potentials=potentials)
        write_potcar(filename=hs, elements=high_spin_elements, potentials=potentials)

        return

    list_low_spin_states  = list(set(ls.name))
    list_high_spin_states = list(set(hs.name))

    if isinstance(ms,str):
        low_spin_length  = len(ls.name)

        if not len(list_low_spin_states) == len(list_high_spin_states): raise RecursionError(f"more than two spin-states detected: ls = {list_low_spin_states}, hs = {list_high_spin_states}")

        ms = permute_spin_staes(ls=list_low_spin_states[0], hs=list_high_spin_states[0], length=low_spin_length, skip_reference=False)

    else:
        ms.append(ls.name)
        ms.append(hs.name)

    for i in ms:
        mixed_spin_state = parent_low_spin_state/i

        if not os.path.exists(mixed_spin_state): os.mkdir(mixed_spin_state)

        if types is not None:
            mixed_spin_magmom = set_magnetic_moments_from_types(spinstate=[i, list_low_spin_states[0]], atoms=low_spin_atoms, elements=low_spin_elements, types=types)
        else:
            mixed_spin_magmom = set_magnetic_moments_from_magmom(spinstate=[i, list_low_spin_states[0]], atoms=low_spin_atoms, elements=low_spin_elements, ls=low_spin_magmom, hs=high_spin_magmom)

        write_incar(filename=mixed_spin_state, magmom=mixed_spin_magmom)
        write_potcar(filename=mixed_spin_state, elements=low_spin_elements, potentials=potentials)

        if not fileio(filename=mixed_spin_state/"POSCAR", option="s"): copy_poscar(source=ls, destiny=mixed_spin_state)

#############################################################################################################################
