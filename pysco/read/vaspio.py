#############################################################################################################################

##       #### ########  ########     ###    ########  #### ########  ######
##        ##  ##     ## ##     ##   ## ##   ##     ##  ##  ##       ##    ##
##        ##  ##     ## ##     ##  ##   ##  ##     ##  ##  ##       ##
##        ##  ########  ########  ##     ## ########   ##  ######    ######
##        ##  ##     ## ##   ##   ######### ##   ##    ##  ##             ##
##        ##  ##     ## ##    ##  ##     ## ##    ##   ##  ##       ##    ##
######## #### ########  ##     ## ##     ## ##     ## #### ########  ######

#############################################################################################################################

import shutil
import textwrap

import numpy as np

from pathlib import Path as path
#
### local libraries
#
from ..thermo.constants import spin_crossover_ions, spin_crossover_metals

from .strings import vasp_total_energy_str, vasp_cell_volume_str, vasp_fermi_energy_str
from .strings import vasp_frequency_str, vasp_mag_start_str

#############################################################################################################################

##     ##  #######  ########  ##     ## ##       ########  ######
###   ### ##     ## ##     ## ##     ## ##       ##       ##    ##
#### #### ##     ## ##     ## ##     ## ##       ##       ##
## ### ## ##     ## ##     ## ##     ## ##       ######    ######
##     ## ##     ## ##     ## ##     ## ##       ##             ##
##     ## ##     ## ##     ## ##     ## ##       ##       ##    ##
##     ##  #######  ########   #######  ######## ########  ######

#############################################################################################################################
### [READ] Get data from OUTCAR

def _readout(data:list=None, string:str="") -> int:

    idx = (idx for idx in np.arange(data.size-1, 0, -1) if string in data[idx])

    return idx

def get_from_outcar(data:list=None)-> tuple[float, float, float]:

    idx_ene = next(_readout(data=data, string=vasp_total_energy_str))
    idx_vol = next(_readout(data=data[:idx_ene], string=vasp_cell_volume_str))
    idx_fer = next(_readout(data=data[:idx_vol], string=vasp_fermi_energy_str))

    total_energy = float(data[idx_ene].split()[6])
    cell_volume  = float(data[idx_vol].split()[4])
    fermi_energy = float(data[idx_fer].split()[2])

    if not total_energy: raise RuntimeError("total energy not found in OUTCAR")
    if not fermi_energy: raise RuntimeError("Fermi energy not found in OUTCAR")

    try:
        cell_volume += 0.0

    except NameError:
        cell_volume = 0.0

    return total_energy, fermi_energy, cell_volume

#############################################################################################################################
### [READ] Get data from EIGENVAL

def get_from_eigcar(data:list=None, fermi:float=0.0) -> np.ndarray:

    details  = data[5].rstrip().split()

    spin     = int(data[0].rstrip().split()[-1])
    kpoints  = int(details[1])
    bands    = int(details[2])

#    points   = [j for i in range(kpoints,len(data), bands+2) for j in range(i,i+bands)]
    points   = [j for i in range(8,len(data), bands+2) for j in range(i,i+bands)]

    alpha    = np.loadtxt(data[points], usecols=1, dtype=float)

    if spin == 2: beta = np.loadtxt(data[points], usecols=2, dtype=float)

    orbital_energy = alpha - fermi if spin == 1 else np.concatenate((alpha,beta)) - fermi
    
    return orbital_energy

#############################################################################################################################
### [READ] Get data from POSCAR

def get_from_poscar(data:list=None) -> tuple[np.ndarray, list]:

    elements = data[5].rstrip().split()
    atoms    = data[6].rstrip().split()

    if any(i.isdigit() for i in elements):  raise RuntimeError("unexpected element symbol in POSCAR")

    if any(not i.isdigit() for i in atoms): raise RuntimeError("unexpected atom number in POSCAR")

    atoms = np.asarray([int(i) for i in atoms])

    return atoms, elements

#############################################################################################################################
### [READ] Get data from VIBCAR

def get_from_vibcar(data:list=None) -> tuple[np.ndarray, float]:

    imaginary = np.char.find(data, sub="f/i")

    if any(i > 0 for i in imaginary): raise RuntimeError("found imaginary frequency in VIBCAR")

    freqdata = np.loadtxt(data, usecols=(3,9), dtype=float)

    frequencies       = freqdata[:,0]*1e12
    
    zero_point_energy = np.sum(freqdata[:,1]/2000)

    return frequencies, zero_point_energy

#############################################################################################################################
### [READ] Get harmonic frequencies from OUTCAR

def _vibindex(data:list=None, vibstr:str=vasp_frequency_str, string:str="Finite differences POTIM"):

    index = []

    for idx, i in enumerate(data):
        if vibstr in i: index.append[idx]
        if string in i: break

    index = np.array(index)

    vibdata = data[np.asarray(index)]

    return np.array(vibdata)

def vib_from_outcar(data:path=None) -> tuple[np.ndarray, float]:

    vibdata = []

    with open(data) as outcar:
        for line in outcar:
            if vibstr in line:
                vibdata.append(line.split())

    vibdata  = np.asarray(vibdata)

    freqdata = np.asarray(vibdata[:, [3,9]], dtype=float)

    #idx = next(_readout(data=data, string="Eigenvectors and eigenvalues of the dynamical matrix"))

    #freqdata = _vibindex(data[idx:])

    ##freqdata = np.loadtxt(vibdata, usecols=(3,9), dtype=float)

    frequencies       = freqdata[:,0]*1e12

    zero_point_energy = np.sum(freqdata[:,1]/2000)

    return frequencies, zero_point_energy

#############################################################################################################################
### [READ] Get magnetization from MAGCAR

def mag_from_magcar(data:list=None, string:str=vasp_mag_start_str) -> np.ndarray:

    for idx in np.arange(data.size-1, 0, -1):
        if string in data[idx]: break
        
    for jdx in np.arange(idx+4, data.size, 1):
        if "----" in data[jdx]: break

    states = len(data[idx+2].split())

    if states < 7: raise RuntimeError("d orbitals not found")

    states = 4 if states == 7 else 5

    magvec = np.loadtxt(data[idx+4:jdx], usecols=states, dtype=float)

    return magvec

#############################################################################################################################
### [READ] Get magnetization from OUTCAR

def mag_from_oszcar(data:list=None, string:str="mag=") -> float:
    
    for idx in np.arange(data.size-1, 0, -1):
        if string in data[idx]: break

    total_mag = float(data[idx].split()[-1])
    
    return total_mag

def mag_from_outcar(data:list=None, magnetization:float=0.0, maxiter:int=10) -> np.ndarray:

    iteration = 0

    magvec = mag_from_magcar(data=data)
    
    if magnetization == 0.0:
        magvec[:] = 0.0
        
        return magvec
    
    for iteration in range(maxiter):
        total = np.sum(magvec)

        if total == magnetization:
            break

        else:
            weight         = magvec/total
            #resize         = np.absolute(weight)<1.0/total
            resize         = np.absolute(weight)<np.sum(weight)/weight.size
            weight[resize] = 0.0
            magvec         = weight*magnetization

    if iteration == maxiter - 1: raise RuntimeError("unable to determine magnetic moments, create MAGCAR instead")

    return magvec

#############################################################################################################################
### [WRITE] Set local magnetic moment for spin-state using types dictionary

def set_magnetic_moments_from_types(spinstate:str=None, atoms:list=None, elements:list=None, types:dict=None) -> np.ndarray:

    idx                    = 0
    local_magnetic_moments = []

    elements_list          = np.repeat(elements, atoms)

    if isinstance(spinstate,list):
        
        mixed_spin_state, low_spin_state = list(spinstate)

        mixed_spin_state = np.array(list(mixed_spin_state))
        mask_spin_state  = mixed_spin_state == low_spin_state

    for element in elements_list:
        magnetic_moment = types.get(element, 0.0)

        if isinstance(magnetic_moment,str): magnetic_moment = spin_crossover_ions.get(magnetic_moment, 0.0)

        if isinstance(magnetic_moment,dict) and not isinstance(spinstate,list):
            magnetic_moment  = magnetic_moment[spinstate]

        if isinstance(magnetic_moment,dict) and isinstance(spinstate,list):
            try:
                which_spin_state = "ls" if mask_spin_state[idx] else "hs"
            except IndexError:
                raise RuntimeError(f"missmatch between number of configurations and atoms")
            
            magnetic_moment  = magnetic_moment[which_spin_state]

            idx += 1

        local_magnetic_moments.append(magnetic_moment)

    return np.asarray(local_magnetic_moments)

#############################################################################################################################
### [WRITE] Set local magnetic moment for spin-state using magmom list

def set_magnetic_moments_from_magmom(spinstate:str=None, ls:list=None, hs:list=None, atoms:list=None, elements:list=None, metals:tuple=spin_crossover_metals) -> np.ndarray:

    idx                    = 0
    local_magnetic_moments = []
    is_spin_crossover      = ls != hs
    elements_list          = np.repeat(elements, atoms)

    spin_crossover_centers = np.count_nonzero(is_spin_crossover)
    low_spin_magenitazion  = np.sum(ls[is_spin_crossover])/spin_crossover_centers
    high_spin_magenitazion = np.sum(hs[is_spin_crossover])/spin_crossover_centers

    mixed_spin_state, low_spin_state = list(spinstate)

    mixed_spin_state = np.array(list(mixed_spin_state))
    mask_spin_state  = mixed_spin_state == low_spin_state
    
    for element in elements_list:
        if element in metals:
            try:
                which_spin_state = low_spin_magenitazion if mask_spin_state[idx] else high_spin_magenitazion
            except IndexError:
                raise RuntimeError(f"missmatch between number of configurations and atoms")
            
            local_magnetic_moments.append(which_spin_state)
            
            idx += 1

        else:
            local_magnetic_moments.append(0.0)

    return local_magnetic_moments

#############################################################################################################################
### [WRITE] Generate a potcar file

def write_potcar(filename:path=None, elements:list=None, potentials:dict=None):

    potentials_path = potentials.get("path", None)

    if potentials_path is None: return
    if isinstance(potentials_path,str): potentials_path = path(potentials_path)

    with open(filename/"POTCAR","wb") as potcar:
        for potential in [potentials.get(i, i) for i in elements]:

            with open(potentials_path/potential/"POTCAR","rb") as pot:
                shutil.copyfileobj(pot, potcar)

#############################################################################################################################
### [WRITE] Generate an incar file 

def write_incar(filename:path=None, magmom:np.ndarray=None) -> None:

    tag = ""
    cnt = 0
    siz = magmom.size

    for idx, i in enumerate(magmom):
        now  = i
        nxt  = None if idx == siz -1 else magmom[idx+1]
        cnt += 1

        if not now == nxt:
            tag = f"{tag} {cnt}*{now}"
            cnt = 0

    incar_tags = {"system":filename.name, "nupdown":np.sum(magmom), "magmom":tag}

    open(filename/"INCAR", "w").write(textwrap.dedent("""\
    SYSTEM   = {system}

    GGA      = PE
    NELM     = 999
    ENCUT    = 600
    EDIFF    = 1.E-08
    LASPH    = .TRUE.
    LREAL    = .FALSE.
    
    ISIF     = 3
    IBRION   = 2
    NSW      = 150
    POTIM    = 0.65
    EDIFFG   = -1.E-03

    ISYM     = 2
    LMAXMIX  = 4
    ALGO     = All
    PREC     = Accurate
    ADDGRID  = .TRUE.

    ISMEAR   = 0
    LORBIT   = 11
    SIGMA    = 1.E-02

    ISPIN    = 2
    NUPDOWN  = {nupdown}
    MAGMOM   ={magmom}

    KSPACING = 0.1
    
    """).format(**incar_tags))

#############################################################################################################################
### [WRITE] Copy poscar file

def copy_poscar(source:path=None, destiny:path=None) -> None:

    shutil.copy(source/"POSCAR", destiny/"POSCAR")

#############################################################################################################################
