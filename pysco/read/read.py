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

import numpy as np

from pathlib import Path as path
#
### local libraries
#
from ..thermo.constants import spin_crossover_metals
from ..thermo.constants import cm1_to_hz, cm1_to_ev, au_to_ev
from ..thermo.constants import h_over_k, c_light, h_bar, bohr_radius, atomic_mass

from .strings import vasp_total_energy_str, vasp_fermi_energy_str
from .strings import vasp_mag_start_str, vasp_mag_stop_str
from .strings import vasp_cell_volume_str
from .strings import vasp_frequency_str, vasp_imag_freq_str

from .strings import gaussian_total_energy_str, gaussian_multip_str
from .strings import gaussian_up_occ_str, gaussian_up_vir_str, gaussian_dw_occ_str, gaussian_dw_vir_str
from .strings import gaussian_pop_start_str, gaussian_pop_stop_str
from .strings import gaussian_frequency_str, gaussian_rot_sym_str, gaussian_rot_temp_str

from .strings import nwchem_single_point_str, nwchem_total_energy_str, nwchem_multip_str
from .strings import nwchem_orb_str, nwchem_up_orb_str, nwchem_dw_orb_str, nwchem_occ_str
from .strings import nwchem_geom_str, nwchem_overlap_str
from .strings import nwchem_frequency_str, nwchem_rot_sym_str, nwchem_rot_temp_str

from .strings import orca_total_energy_str, orca_multip_str, orca_general_stop_str
from .strings import orca_orbital_str, orca_up_orb_str, orca_dw_orb_str
from .strings import orca_pop_start_str, orca_pop_stop_str
from .strings import orca_frequency_str, orca_rot_sym_str, orca_rot_const_str

#############################################################################################################################
### ATTENTION: Information regarding units from output files

# energy          = > eV
# volume          = > Angstroem^3
# jahnteller      = > dimensionless
# frequencies     = > Hz
# magnetization   = > up-electrons - down-electrons
# orbitalenergies = > eV
# zeropointenergy = > eV

#############################################################################################################################

##     ##    ###     ######  ########  
##     ##   ## ##   ##    ## ##     ## 
##     ##  ##   ##  ##       ##     ## 
##     ## ##     ##  ######  ########  
 ##   ##  #########       ## ##        
  ## ##   ##     ## ##    ## ##        
   ###    ##     ##  ######  ##

#############################################################################################################################

class vasp:

    __slots__ = (
                 ### MANDATORY attributes:
                 "atoms", "orbit", "energy", "volume", "elements", "rotationalsymm", "rotationaltemp",
                 "jahnteller", "frequencies", "magnetization", "orbitalenergies", "zeropointenergy", "step"
                )
#
### Initialize
#
    def __init__(self, filesdir:str|path,
                 magnetization:list=None, jahnteller:float=1.0, orbit:float=1.0) -> None:
        
        ### start MANDATORY attributes
        self.atoms           = None
        self.orbit           = None
        self.energy          = None
        self.volume          = None
        self.elements        = None
        self.jahnteller      = None
        self.frequencies     = None
        self.magnetization   = None
        self.orbitalenergies = None
        self.zeropointenergy = None
        
        self.rotationalsymm  = None
        self.rotationaltemp  = None

        self.step            = []
        ### end MANDATORY attributes
        
        if isinstance(filesdir,str): filesdir = path(filesdir)

        try:
            self.jahnteller = float(jahnteller)
        except:
            raise RuntimeError(f"jahnteller = {jahnteller} is not a number")
        
        try:
            self.orbit = float(orbit)
        except:
            raise RuntimeError(f"orbit = {orbit} is not a number")
        
        if (magnetization is not None) and isinstance(magnetization, (list, np.ndarray)):
            self.magnetization = np.asarray(magnetization)
#
### MANDATORY: Extract [atoms, energy, volume, elements, frequencies, magnetization, orbitalenergies, zeropointenergy]
#
        magnetization      = self._get_output_data(filesdir)
        self.magnetization = self._set_magnetic_moments(magnetization)
#
### MANDATORY: Prepare tuple for use as args
#
    def zipit(self, interaction:float=0.0, nhs:float=None) -> tuple:

        if not self.step: raise RuntimeError("spin-step is empty")
        
        crossover_energy, spin_entropy, volume_expansion, frequencies, orbitalenergies, rotational = self.step

        nhs_tag    = ("nhs",         nhs              )
        sco_tag    = ("E_sco",       crossover_energy )
        s_spin_tag = ("S_spin",      spin_entropy     )
        pdv_tag    = ("PdV",         volume_expansion )
        freq_tag   = ("frq",         frequencies      )
        eig_tag    = ("eig",         orbitalenergies  )
        inter_tag  = ("interaction", interaction      )
        rot_tag    = ("rot",         rotational       )
        
        if nhs is not None:
            return (nhs_tag, sco_tag, s_spin_tag, pdv_tag, freq_tag, eig_tag, rot_tag, inter_tag)
        
        else:
            return dict(i for i in (sco_tag, s_spin_tag, pdv_tag, freq_tag, rot_tag, eig_tag))
#
### Extract all data from output directory
#
    def _get_output_data(self, filesdir) -> list:
        
        outfile = filesdir/"OUTCAR"
        posfile = filesdir/"POSCAR"
        eigfile = filesdir/"EIGENVAL"

        magfile = filesdir/"MAGCAR"
        vibfile = filesdir/"VIBCAR"
        
        fileio(outfile, option="r")
        fileio(posfile, option="r")
        fileio(eigfile, option="r")

        skipped_lines       = 0

        magnetization       = []
        frequencies         = []
        #
        ### Runtime variables
        #
        found_magnetization = False
        skip_mag_lines      = 3
        #
        ### Total Energy, Fermi Energy, and Volume
        #
        with open(outfile, "r") as out:
            for line in out:
                line = line.strip()

                if vasp_cell_volume_str  in line: self.volume  = float(line.split()[4])
                if vasp_fermi_energy_str in line: fermi_energy = float(line.split()[2])
                if vasp_total_energy_str in line: self.energy  = float(line.split()[4])
        #
        ### Elements
        #
        with open(posfile, "r") as out:
            for idx, line in enumerate(out):
                if idx == 5:
                    self.elements = np.asarray(line.strip().split(), dtype=str)
                    continue
                
                if idx == 6:
                    self.atoms    = np.asarray(line.strip().split(), dtype=int)
                    break
        #
        ### Orbital Energies
        #
        data = np.asarray( open(eigfile).readlines() )

        details  = data[5].rstrip().split()
        spin     = int(data[0].rstrip().split()[-1])
        kpoints  = int(details[1])
        bands    = int(details[2])

        points   = [j for i in range(8,len(data), bands+2) for j in range(i,i+bands)]

        alpha    = np.loadtxt(data[points], usecols=1, dtype=float)

        if spin == 2: beta = np.loadtxt(data[points], usecols=2, dtype=float)

        self.orbitalenergies = alpha - fermi_energy if spin == 1 else np.concatenate((alpha,beta)) - fermi_energy
        #
        ### Harmonic Frequencies
        #
        if fileio(vibfile, option="s"):
            with open(vibfile, "r") as out:
                for line in out:
                    line = line.strip()

                    if vasp_frequency_str in line:
                        if vasp_imag_freq_str in line: raise RuntimeError(f"found imaginary frequency in {vibfile.name}")
                        frequencies += [float(line.split()[7])]
        #
        ### Magnetization
        #
        if (self.magnetization is None) and fileio(magfile, option="s"):
            with open(magfile, "r") as out:
                for line in out:
                    line = line.strip()

                    if vasp_mag_start_str in line:
                        skipped_lines       = 0
                        found_magnetization = True
                        continue

                    if   (found_magnetization) and (skipped_lines < skip_mag_lines):
                        skipped_lines += 1
                        continue
                    elif (found_magnetization) and (vasp_mag_stop_str in line):
                        skipped_lines = 0
                        found_magnetization = False
                        continue
                    elif (found_magnetization) and (skipped_lines == skip_mag_lines):
                        magnetization += [ float(line.split()[-1]) ]
                        continue
        
        if not frequencies: frequencies = 3*[0.0]
        
        frequencies = np.asarray(frequencies, dtype=float)
        frequencies = frequencies[frequencies != 0.0]
        
        self.frequencies     = frequencies*cm1_to_hz
        self.zeropointenergy = 0.5*np.sum(frequencies)*cm1_to_ev

        return magnetization
#
### Set magnetic moments
#
    def _set_magnetic_moments(self, magnetization):

        if self.magnetization is not None: return self.magnetization
    
        if not magnetization: raise RuntimeError("either set 'magnetization=' or provide a MAGCAR")

        return np.asarray(magnetization, dtype=float)
        
#############################################################################################################################

########  ######  ########  ########  ########  ######   ######   #######  
##       ##    ## ##     ## ##     ## ##       ##    ## ##    ## ##     ## 
##       ##       ##     ## ##     ## ##       ##       ##       ##     ## 
######    ######  ########  ########  ######    ######   ######  ##     ## 
##             ## ##        ##   ##   ##             ##       ## ##     ## 
##       ##    ## ##        ##    ##  ##       ##    ## ##    ## ##     ## 
########  ######  ##        ##     ## ########  ######   ######   #######

#############################################################################################################################


#############################################################################################################################

 ######      ###    ##     ##  ######   ######  ####    ###    ##    ## 
##    ##    ## ##   ##     ## ##    ## ##    ##  ##    ## ##   ###   ## 
##         ##   ##  ##     ## ##       ##        ##   ##   ##  ####  ## 
##   #### ##     ## ##     ##  ######   ######   ##  ##     ## ## ## ## 
##    ##  ######### ##     ##       ##       ##  ##  ######### ##  #### 
##    ##  ##     ## ##     ## ##    ## ##    ##  ##  ##     ## ##   ### 
 ######   ##     ##  #######   ######   ######  #### ##     ## ##    ##

#############################################################################################################################

class gaussian:

    __slots__ = (
                 ### MANDATORY attributes:
                 "atoms", "orbit", "energy", "volume", "elements", "rotationalsymm", "rotationaltemp",
                 "jahnteller", "frequencies", "magnetization", "orbitalenergies", "zeropointenergy", "step"
                )
#
### Initialize
#
    def __init__(self, outfile:str|path,
                 magnetization:list=None, jahnteller:float=1.0, orbit:float=1.0) -> None:
        
        ### start MANDATORY attributes
        self.atoms           = None
        self.orbit           = None
        self.energy          = None
        self.volume          = None
        self.elements        = None
        self.jahnteller      = None
        self.frequencies     = None
        self.magnetization   = None
        self.orbitalenergies = None
        self.zeropointenergy = None
        
        self.rotationalsymm  = None
        self.rotationaltemp  = None

        self.step            = []
        ### end MANDATORY attributes
        
        if isinstance(outfile,str): outfile = path(outfile)
        fileio(outfile, option="r")

        try:
            self.jahnteller = float(jahnteller)
        except:
            raise RuntimeError(f"jahnteller = {jahnteller} is not a number")
        
        try:
            self.orbit = float(orbit)
        except:
            raise RuntimeError(f"orbit = {orbit} is not a number")
        
        if (magnetization is not None) and isinstance(magnetization, (list, np.ndarray)):
            self.magnetization = np.asarray(magnetization)
#
### MANDATORY: Extract [atoms, energy, volume, elements, frequencies, magnetization, orbitalenergies, zeropointenergy]
#
        multiplicity       = self._get_output_data(outfile)
        self.magnetization = self._set_magnetic_moments(multiplicity)
#
### MANDATORY: Prepare tuple for use as args
#
    def zipit(self, interaction:float=0.0, nhs:float=None) -> tuple:

        if not self.step: raise RuntimeError("spin-step is empty")
        
        crossover_energy, spin_entropy, volume_expansion, frequencies, orbitalenergies, rotational = self.step

        nhs_tag    = ("nhs",         nhs              )
        sco_tag    = ("E_sco",       crossover_energy )
        s_spin_tag = ("S_spin",      spin_entropy     )
        pdv_tag    = ("PdV",         volume_expansion )
        freq_tag   = ("frq",         frequencies      )
        eig_tag    = ("eig",         orbitalenergies  )
        inter_tag  = ("interaction", interaction      )
        rot_tag    = ("rot",         rotational       )
        
        if nhs is not None:
            return (nhs_tag, sco_tag, s_spin_tag, pdv_tag, freq_tag, eig_tag, rot_tag, inter_tag)
        
        else:
            return dict(i for i in (sco_tag, s_spin_tag, pdv_tag, freq_tag, eig_tag, rot_tag))
#
### Extract all data from output file
#
    def _get_output_data(self, data) -> float:
        
        skipped_lines        = 0

        multiplicity         = None

        elements             = []
        frequencies          = []
        #
        ### Runtime variables
        #
        found_population     = False
        skip_pop_lines       = 1

        with open(data, "r") as out:
            for line in out:
                line = line.strip()
                #
                ### Total Energy and Multiplicity
                #
                if gaussian_multip_str in line: multiplicity = float(line.split()[5])
                if gaussian_total_energy_str in line:
                    self.energy = float(line.split()[4])*au_to_ev
                    
                    alpha_orbital_energy = []
                    beta_orbital_energy  = []
                    alpha_orbital_occup  = []
                    beta_orbital_occup   = []
                #
                ### Orbital Energies and Occupations
                #
                if (gaussian_up_occ_str in line) or (gaussian_dw_occ_str in line):
                    energy     = [ float(i) for i in line.split()[4:] ]
                    occupation = [ 1.0 for i in line.split()[4:] ]

                    if gaussian_up_occ_str in line:
                        alpha_orbital_energy += energy
                        alpha_orbital_occup  += occupation
                    elif gaussian_dw_occ_str in line:
                        beta_orbital_energy  += energy
                        beta_orbital_occup   += occupation

                elif (gaussian_up_vir_str in line) or (gaussian_dw_vir_str in line):
                    energy     = [ float(i) for i in line.split()[4:] ]
                    occupation = [ 0.0 for i in line.split()[4:] ]

                    if gaussian_up_vir_str in line:
                        alpha_orbital_energy += energy
                        alpha_orbital_occup  += occupation
                    elif gaussian_dw_vir_str in line:
                        beta_orbital_energy  += energy
                        beta_orbital_occup   += occupation
                #
                ### Elements
                #
                if gaussian_pop_start_str in line:
                    skipped_lines    = 0
                    found_population = True
                    continue

                if   (found_population) and (skipped_lines < skip_pop_lines):
                    skipped_lines   += 1
                    continue
                elif (found_population) and (gaussian_pop_stop_str in line):
                    skipped_lines    = 0
                    found_population = False
                    continue
                elif (found_population) and (skipped_lines == skip_pop_lines):
                    elements += [ line.split()[1] ]
                    continue
                #
                ### Harmonic Frequencies
                #
                if gaussian_frequency_str in line: frequencies += [ float(i) for i in line.split()[2:] ]
                #
                ### Rotational data
                #
                if gaussian_rot_sym_str  in line:
                    self.rotationalsymm = 1.0/float(line.split()[3])
                    
                if gaussian_rot_temp_str in line:
                    self.rotationaltemp = np.pi/np.prod([ float(i) for i in line.split()[3:] ])
        
        alpha_orbital_energy = np.asarray(alpha_orbital_energy, dtype=float)*au_to_ev
        beta_orbital_energy  = np.asarray(beta_orbital_energy,  dtype=float)*au_to_ev
        alpha_orbital_occup  = np.asarray(alpha_orbital_occup,  dtype=float)
        beta_orbital_occup   = np.asarray(beta_orbital_occup,   dtype=float)

        if beta_orbital_occup.size:
            fermi_energy = max(alpha_orbital_energy[alpha_orbital_occup==1.0][-1], beta_orbital_energy[beta_orbital_occup==1.0][-1])
            self.orbitalenergies = np.concatenate((alpha_orbital_energy, beta_orbital_energy)) - fermi_energy

        else:
            fermi_energy = alpha_orbital_energy[alpha_orbital_occup==1.0][-1]
            self.orbitalenergies = np.concatenate((alpha_orbital_energy, alpha_orbital_energy)) - fermi_energy
    
        self.elements = np.asarray(elements,   dtype=str)
        self.atoms    = np.ones_like(elements, dtype=int)
        
        if not frequencies: frequencies = 3*[0.0]
        if any(i < 0 for i in frequencies): raise RuntimeError(f"found imaginary frequency in {data.name}")

        frequencies = np.asarray(frequencies, dtype=float)
        frequencies = frequencies[frequencies != 0.0]
        
        self.frequencies     = frequencies*cm1_to_hz
        self.zeropointenergy = 0.5*np.sum(frequencies)*cm1_to_ev
        
        self.volume = 0.0

        return multiplicity
#
### Set magnetic moments per atom
#
    def _set_magnetic_moments(self, multiplicity):

        if self.magnetization is not None:
            return self.magnetization

        if multiplicity == 0.0:
            raise RuntimeError(f"cannot resolve magnetic moments. Use 'magnetization=' instead")

        magnetization = np.array([multiplicity-1.0 if i in spin_crossover_metals else 0.0 for i in self.elements], dtype=float)

        return magnetization

#############################################################################################################################

##    ## ##      ##  ######  ##     ## ######## ##     ## 
###   ## ##  ##  ## ##    ## ##     ## ##       ###   ### 
####  ## ##  ##  ## ##       ##     ## ##       #### #### 
## ## ## ##  ##  ## ##       ######### ######   ## ### ## 
##  #### ##  ##  ## ##       ##     ## ##       ##     ## 
##   ### ##  ##  ## ##    ## ##     ## ##       ##     ## 
##    ##  ###  ###   ######  ##     ## ######## ##     ##

#############################################################################################################################

class nwchem:

    __slots__ = (
                 ### MANDATORY attributes:
                 "atoms", "orbit", "energy", "volume", "elements", "rotationalsymm", "rotationaltemp",
                 "jahnteller", "frequencies", "magnetization", "orbitalenergies", "zeropointenergy", "step"
                )
#
### Initialize
#
    def __init__(self, outfile:str|path,
                 magnetization:list=None, jahnteller:float=1.0, orbit:float=1.0) -> None:
        
        ### start MANDATORY attributes
        self.atoms           = None
        self.orbit           = None
        self.energy          = None
        self.volume          = None
        self.elements        = None
        self.jahnteller      = None
        self.frequencies     = None
        self.magnetization   = None
        self.orbitalenergies = None
        self.zeropointenergy = None
        
        self.rotationalsymm  = None
        self.rotationaltemp  = None

        self.step            = []
        ### end MANDATORY attributes
        
        if isinstance(outfile,str): outfile = path(outfile)
        fileio(outfile, option="r")

        try:
            self.jahnteller = float(jahnteller)
        except:
            raise RuntimeError(f"jahnteller = {jahnteller} is not a number")
        
        try:
            self.orbit = float(orbit)
        except:
            raise RuntimeError(f"orbit = {orbit} is not a number")
        
        if (magnetization is not None) and isinstance(magnetization, (list, np.ndarray)):
            self.magnetization = np.asarray(magnetization)
#
### MANDATORY: Extract [atoms, energy, volume, elements, frequencies, magnetization, orbitalenergies, zeropointenergy]
#
        multiplicity       = self._get_output_data(outfile)
        self.magnetization = self._set_magnetic_moments(multiplicity)        
#
### MANDATORY: Prepare tuple for use as args
#
    def zipit(self, interaction:float=0.0, nhs:float=None) -> tuple:

        if not self.step: raise RuntimeError("spin-step is empty")
        
        crossover_energy, spin_entropy, volume_expansion, frequencies, orbitalenergies, rotational = self.step

        nhs_tag    = ("nhs",         nhs              )
        sco_tag    = ("E_sco",       crossover_energy )
        s_spin_tag = ("S_spin",      spin_entropy     )
        pdv_tag    = ("PdV",         volume_expansion )
        freq_tag   = ("frq",         frequencies      )
        eig_tag    = ("eig",         orbitalenergies  )
        inter_tag  = ("interaction", interaction      )
        rot_tag    = ("rot",         rotational       )
        
        if nhs is not None:
            return (nhs_tag, sco_tag, s_spin_tag, pdv_tag, freq_tag, eig_tag, rot_tag, inter_tag)
        
        else:
            return dict(i for i in (sco_tag, s_spin_tag, pdv_tag, freq_tag, eig_tag, rot_tag))
#
### Extract all data from output file
#
    def _get_output_data(self, data) -> float:
        
        skipped_lines           = 0

        multiplicity            = None
        total_energy            = None

        alpha_orbital_energy    = []
        beta_orbital_energy     = []
        alpha_orbital_occup     = []
        beta_orbital_occup      = []

        elements                = []
        frequencies             = []
        rotational_temperatures = []
        #
        ### Runtime variables
        #
        found_single_point      = False
        found_spin_up           = False
        found_spin_dw           = False

        found_geometry          = False
        skip_geom_lines         = 3

        found_rotational        = False
        skip_rotational_lines   = 1

        with open(data, "r") as out:
            for line in out:
                line = line.strip()
                #
                ### Elements
                #
                if nwchem_geom_str in line:
                    skipped_lines  = 0
                    found_geometry = True
                    continue

                if   (found_geometry) and (skipped_lines < skip_geom_lines):
                    skipped_lines += 1
                    continue

                if   (found_geometry) and (not line):
                    skipped_lines  = 0
                    found_geometry = False
                    continue
                    
                elif (found_geometry) and (skipped_lines == skip_geom_lines):
                    elements += [line.split()[0]]
                #
                ### Total Energy and Multiplicity
                #
                if nwchem_single_point_str in line:
                    found_single_point = True
                    continue

                if (nwchem_multip_str       in line) and (multiplicity is None): multiplicity = float(line.split()[2])
                if (nwchem_total_energy_str in line) and (total_energy is None): self.energy = float(line.split()[4])*au_to_ev
                #
                ### Orbital Energies and Occupations
                #
                if (found_single_point) and ( (nwchem_up_orb_str in line) or (nwchem_orb_str in line) ):
                    found_spin_up = True
                    found_spin_dw = False
                    continue

                if (found_single_point) and (nwchem_dw_orb_str in line):
                    found_spin_up = False
                    found_spin_dw = True
                    continue

                if (found_single_point) and (nwchem_overlap_str in line):
                    found_spin_up      = False
                    found_spin_dw      = False
                    found_single_point = False
                    continue

                if (found_single_point) and (nwchem_occ_str in line):
                    data_line       = line.split()
                    negative_energy = len(data_line) < 5

                    occupation_data, energy_data = data_line[2:] if negative_energy else data_line[2::2]

                    energy     = float(energy_data.partition("=")[2].replace("D","e")) if negative_energy else float(energy_data.replace("D","e"))
                    occupation = float(occupation_data.partition("=")[2].replace("D","e"))

                    if   found_spin_up:
                        alpha_orbital_energy += [energy]
                        alpha_orbital_occup  += [occupation]
                        continue
                        
                    elif found_spin_dw:
                        beta_orbital_energy  += [energy]
                        beta_orbital_occup   += [occupation]
                        continue
                #
                ### Rotational data
                #
                if nwchem_rot_sym_str in line:
                    self.rotationalsymm = 1.0/float(line.split()[8].replace(")",""))
                #    
                if nwchem_rot_temp_str in line:
                    skipped_lines    = 0
                    found_rotational = True
                    continue

                if   (found_rotational) and (skipped_lines < skip_rotational_lines):
                    skipped_lines += 1
                    continue
                    
                elif (found_rotational) and (not line):
                    skipped_lines       = 0
                    found_rotational    = False
                    self.rotationaltemp = np.pi/np.prod(rotational_temperatures)
                    continue
                    
                elif (found_rotational) and (skipped_lines == skip_rotational_lines):
                    rotational_temperatures += [ float(line.split()[4]) ]
                #
                ### Harmonic Frequencies
                #
                if nwchem_frequency_str in line:
                    frequencies += [ float(i) for i in line.split()[1:] ]

        alpha_orbital_energy = np.asarray(alpha_orbital_energy, dtype=float)*au_to_ev
        beta_orbital_energy  = np.asarray(beta_orbital_energy,  dtype=float)*au_to_ev
        alpha_orbital_occup  = np.asarray(alpha_orbital_occup,  dtype=float)
        beta_orbital_occup   = np.asarray(beta_orbital_occup,   dtype=float)

        if beta_orbital_occup.size:
            fermi_energy = max(alpha_orbital_energy[alpha_orbital_occup==1.0][-1], beta_orbital_energy[beta_orbital_occup==1.0][-1])
            self.orbitalenergies = np.concatenate((alpha_orbital_energy, beta_orbital_energy)) - fermi_energy

        else:
            fermi_energy = alpha_orbital_energy[alpha_orbital_occup==1.0][-1]
            self.orbitalenergies = np.concatenate((alpha_orbital_energy, alpha_orbital_energy)) - fermi_energy
    
        self.elements = np.asarray(elements,   dtype=str)
        self.atoms    = np.ones_like(elements, dtype=int)
        
        if not frequencies: frequencies = 3*[0.0]
        if any(i < 0 for i in frequencies): raise RuntimeError(f"found imaginary frequency in {data.name}")

        frequencies = np.asarray(frequencies, dtype=float)
        frequencies = frequencies[frequencies != 0.0]
        
        self.frequencies     = frequencies*cm1_to_hz
        self.zeropointenergy = 0.5*np.sum(frequencies)*cm1_to_ev
        
        self.volume = 0.0
    
        return multiplicity
#
### Set magnetic moments per atom
#
    def _set_magnetic_moments(self, multiplicity):

        if self.magnetization is not None:
            return self.magnetization

        if multiplicity == 0.0:
            raise RuntimeError(f"cannot resolve magnetic moments. Use 'magnetization=' instead")

        magnetization = np.array([multiplicity-1.0 if i in spin_crossover_metals else 0.0 for i in self.elements], dtype=float)

        return magnetization

#############################################################################################################################

 #######  ########   ######     ###    
##     ## ##     ## ##    ##   ## ##   
##     ## ##     ## ##        ##   ##  
##     ## ########  ##       ##     ## 
##     ## ##   ##   ##       ######### 
##     ## ##    ##  ##    ## ##     ## 
 #######  ##     ##  ######  ##     ##

#############################################################################################################################

class orca:

    __slots__ = (
                 ### MANDATORY attributes:
                 "atoms", "orbit", "energy", "volume", "elements", "rotationalsymm", "rotationaltemp",
                 "jahnteller", "frequencies", "magnetization", "orbitalenergies", "zeropointenergy", "step"
                )
#
### Initialize
#
    def __init__(self, outfile:str|path,
                 magnetization:list=None, jahnteller:float=1.0, orbit:float=1.0) -> None:
        
        ### start MANDATORY attributes
        self.atoms           = None
        self.orbit           = None
        self.energy          = None
        self.volume          = None
        self.elements        = None
        self.jahnteller      = None
        self.frequencies     = None
        self.magnetization   = None
        self.orbitalenergies = None
        self.zeropointenergy = None
        
        self.rotationalsymm  = None
        self.rotationaltemp  = None

        self.step            = []
        ### end MANDATORY attributes
        
        if isinstance(outfile,str): outfile = path(outfile)
        fileio(outfile, option="r")

        try:
            self.jahnteller = float(jahnteller)
        except:
            raise RuntimeError(f"jahnteller = {jahnteller} is not a number")
        
        try:
            self.orbit = float(orbit)
        except:
            raise RuntimeError(f"orbit = {orbit} is not a number")
        
        if (magnetization is not None) and isinstance(magnetization, (list, np.ndarray)):
            self.magnetization = np.asarray(magnetization)
#
### MANDATORY: Extract [atoms, energy, volume, elements, frequencies, magnetization, orbitalenergies, zeropointenergy]
#
        multiplicity       = self._get_output_data(outfile)
        self.magnetization = self._set_magnetic_moments(multiplicity)
#
### MANDATORY: Prepare tuple for use as args
#
    def zipit(self, interaction:float=0.0, nhs:float=None) -> tuple:

        if not self.step: raise RuntimeError("spin-step is empty")
        
        crossover_energy, spin_entropy, volume_expansion, frequencies, orbitalenergies, rotational = self.step

        nhs_tag    = ("nhs",         nhs              )
        sco_tag    = ("E_sco",       crossover_energy )
        s_spin_tag = ("S_spin",      spin_entropy     )
        pdv_tag    = ("PdV",         volume_expansion )
        freq_tag   = ("frq",         frequencies      )
        eig_tag    = ("eig",         orbitalenergies  )
        inter_tag  = ("interaction", interaction      )
        rot_tag    = ("rot",         rotational       )
        
        if nhs is not None:
            return (nhs_tag, sco_tag, s_spin_tag, pdv_tag, freq_tag, eig_tag, rot_tag, inter_tag)
        
        else:
            return dict(i for i in (sco_tag, s_spin_tag, pdv_tag, freq_tag, eig_tag, rot_tag))
#
### Extract all data from output file
#
    def _get_output_data(self, data) -> float:

        skipped_lines          = 0

        alpha_orbital_energy   = []
        beta_orbital_energy    = []
        alpha_orbital_occup    = []
        beta_orbital_occup     = []

        elements               = []
        frequencies            = []
        #
        ### Runtime variables
        #
        found_orbital_energies = False
        found_spin_up          = False
        found_spin_dw          = False
        skip_orbital_lines     = 2

        found_population       = False
        skip_pop_lines         = 1

        found_frequencies      = False
        skip_freq_lines        = 4

        with open(data, "r") as out:
            for line in out:
                line = line.strip()
                #
                ### Total Energy and Multiplicity
                #
                if orca_multip_str       in line: multiplicity = float(line.split()[2])
                if orca_total_energy_str in line: self.energy  = float(line.split()[4])*au_to_ev
                #
                ### Orbital Energies and Occupations
                #
                if orca_orbital_str in line:
                    skipped_lines          = 0
                    found_orbital_energies = True
                    continue

                if   (orca_up_orb_str in line) and (found_orbital_energies):
                    skipped_lines += 1
                    found_spin_up  = True
                    continue
                elif (orca_dw_orb_str in line) and (found_orbital_energies):
                    skipped_lines += 1
                    found_spin_dw  = True
                    continue

                elif (found_spin_up) and (found_orbital_energies) and (skipped_lines < skip_orbital_lines):
                    skipped_lines += 1
                    continue
                elif (found_spin_dw) and (found_orbital_energies) and (skipped_lines < skip_orbital_lines):
                    skipped_lines += 1
                    continue

                elif (found_spin_up) and (found_orbital_energies) and ( (not line) or (orca_general_stop_str in line) ):
                    skipped_lines = 0
                    found_spin_up = False
                    continue
                elif (found_spin_dw) and (found_orbital_energies) and ( (not line) or (orca_general_stop_str in line) ):
                    skipped_lines = 0
                    found_spin_dw = False
                    continue

                elif (found_spin_up) and (found_orbital_energies) and (skipped_lines == skip_orbital_lines):
                    alpha_orbital_energy += [ float(line.split()[3])]
                    alpha_orbital_occup  += [ float(line.split()[1])]
                elif (found_spin_dw) and (found_orbital_energies) and (skipped_lines == skip_orbital_lines):
                    beta_orbital_energy  += [ float(line.split()[3])]
                    beta_orbital_occup   += [ float(line.split()[1])]
                #
                ### Elements
                #
                if orca_pop_start_str in line:
                    skipped_lines    = 0
                    found_population = True
                    continue

                if   (found_population) and (skipped_lines < skip_pop_lines):
                    skipped_lines += 1
                    continue
                elif (found_population) and (orca_pop_stop_str in line):
                    skipped_lines    = 0
                    found_population = False
                    continue
                elif (found_population) and (skipped_lines == skip_pop_lines):
                    elements += [ line.split()[1].replace(":","") ]
                #
                ### Harmonic Frequencies
                #
                if orca_frequency_str in line:
                    skipped_lines     = 0
                    found_frequencies = True
                    continue

                if   (found_frequencies) and (skipped_lines < skip_freq_lines):
                    skipped_lines += 1
                    continue
                elif (found_frequencies) and ( (not line) or (orca_general_stop_str in line) ):
                    skipped_lines     = 0
                    found_frequencies = False
                    continue
                elif (found_frequencies) and (skipped_lines == skip_freq_lines):
                    frequencies += [ float(line.split()[1]) ]
                #
                ### Rotational data
                #
                if orca_rot_sym_str in line:
                    self.rotationalsymm = 1.0/float(line.split()[5])

                if orca_rot_const_str in line:
                    rotational_temperatures = np.asarray(line.split()[4:], dtype=float)*1e+02*h_over_k*c_light
                    self.rotationaltemp     = np.pi/np.prod(rotational_temperatures)
        
        alpha_orbital_energy = np.asarray(alpha_orbital_energy, dtype=float)
        beta_orbital_energy  = np.asarray(beta_orbital_energy,  dtype=float)
        alpha_orbital_occup  = np.asarray(alpha_orbital_occup,  dtype=float)
        beta_orbital_occup   = np.asarray(beta_orbital_occup,   dtype=float)

        if beta_orbital_occup.size:
            fermi_energy = max(alpha_orbital_energy[alpha_orbital_occup==1.0][-1], beta_orbital_energy[beta_orbital_occup==1.0][-1])
            self.orbitalenergies = np.concatenate((alpha_orbital_energy, beta_orbital_energy)) - fermi_energy

        else:
            fermi_energy = alpha_orbital_energy[alpha_orbital_occup==1.0][-1]
            self.orbitalenergies = np.concatenate((alpha_orbital_energy, alpha_orbital_energy)) - fermi_energy
    
        self.elements = np.asarray(elements,   dtype=str)
        self.atoms    = np.ones_like(elements, dtype=int)
        
        if not frequencies: frequencies = 3*[0.0]
        if any(i < 0 for i in frequencies): raise RuntimeError(f"found imaginary frequency in {data.name}")

        frequencies = np.asarray(frequencies, dtype=float)
        frequencies = frequencies[frequencies != 0.0]
        
        self.frequencies     = frequencies*cm1_to_hz
        self.zeropointenergy = 0.5*np.sum(frequencies)*cm1_to_ev
        
        self.volume = 0.0
    
        return multiplicity
#
### Set magnetic moments per atom
#
    def _set_magnetic_moments(self, multiplicity):

        if self.magnetization is not None: return self.magnetization

        if multiplicity <= 0.0:
            raise RuntimeError(f"cannot resolve magnetic moments. Use 'magnetization=' instead")

        magnetization = np.array([multiplicity-1.0 if i in spin_crossover_metals else 0.0 for i in self.elements], dtype=float)

        return magnetization

#############################################################################################################################

########  ##    ##  ######   ######  ########
##     ##  ##  ##  ##    ## ##    ## ##
##     ##   ####   ##       ##       ##
########     ##     ######  ##       ######
##           ##          ## ##       ##
##           ##    ##    ## ##    ## ##
##           ##     ######   ######  ##

#############################################################################################################################

class pyscf:

    __slots__ = (
                 ### MANDATORY attributes:
                 "atoms", "orbit", "energy", "volume", "elements", "rotationalsymm", "rotationaltemp",
                 "jahnteller", "frequencies", "magnetization", "orbitalenergies", "zeropointenergy", "step"
                )
#
### Initialize
#
    def __init__(self, meanfield, harmonic_analysis:dict=None,
                 magnetization:list=None, jahnteller:float=1.0, orbit:float=1.0) -> None:
        
        ### start MANDATORY attributes
        self.atoms           = None
        self.orbit           = None
        self.energy          = None
        self.volume          = None
        self.elements        = None
        self.jahnteller      = None
        self.frequencies     = None
        self.magnetization   = None
        self.orbitalenergies = None
        self.zeropointenergy = None
        
        self.rotationalsymm  = None
        self.rotationaltemp  = None

        self.step            = []
        ### end MANDATORY attributes

        try:
            self.jahnteller = float(jahnteller)
        except:
            raise RuntimeError(f"jahnteller = {jahnteller} is not a number")
        
        try:
            self.orbit = float(orbit)
        except:
            raise RuntimeError(f"orbit = {orbit} is not a number")
        
        if (magnetization is not None) and isinstance(magnetization, (list, np.ndarray)):
            self.magnetization = np.asarray(magnetization)
#
### MANDATORY: Extract [atoms, energy, volume, elements, frequencies, magnetization, orbitalenergies, zeropointenergy]
#
        multiplicity       = self._get_output_data(meanfield, harmonic_analysis)
        self.magnetization = self._set_magnetic_moments(multiplicity)
#
### MANDATORY: Prepare tuple for use as args
#
    def zipit(self, interaction:float=0.0, nhs:float=None) -> tuple:

        if not self.step: raise RuntimeError("spin-step is empty")
        
        crossover_energy, spin_entropy, volume_expansion, frequencies, orbitalenergies, rotational = self.step

        nhs_tag    = ("nhs",         nhs              )
        sco_tag    = ("E_sco",       crossover_energy )
        s_spin_tag = ("S_spin",      spin_entropy     )
        pdv_tag    = ("PdV",         volume_expansion )
        freq_tag   = ("frq",         frequencies      )
        eig_tag    = ("eig",         orbitalenergies  )
        inter_tag  = ("interaction", interaction      )
        rot_tag    = ("rot",         rotational       )
        
        if nhs is not None:
            return (nhs_tag, sco_tag, s_spin_tag, pdv_tag, freq_tag, eig_tag, rot_tag, inter_tag)
        
        else:
            return dict(i for i in (sco_tag, s_spin_tag, pdv_tag, freq_tag, eig_tag, rot_tag))
#
### Extract all data from output file
#
    def _get_output_data(self, meanfield, harmonic) -> float:

        alpha_orbital_energy   = []
        beta_orbital_energy    = []
        alpha_orbital_occup    = []
        beta_orbital_occup     = []

        frequencies            = []
        #
        ### Total Energy and Multiplicity
        #
        self.energy  = float(meanfield.e_tot)
        multiplicity = meanfield.mol.multiplicity
        #
        ### Orbital Energies and Occupations
        #
        if meanfield.mo_occ.ndim == 1:
            alpha_orbital_energy = np.asarray(meanfield.mo_energy, dtype=float)
            alpha_orbital_occup  = np.asarray(meanfield.mo_occ, dtype=float)
        else:
            alpha_orbital_energy = np.asarray(meanfield.mo_energy[0], dtype=float)
            beta_orbital_energy  = np.asarray(meanfield.mo_energy[1], dtype=float)
            alpha_orbital_occup  = np.asarray(meanfield.mo_occ[0], dtype=float)
            beta_orbital_occup   = np.asarray(meanfield.mo_occ[1], dtype=float)
        #
        ### Elements
        #
        self.elements = np.asarray(meanfield.mol.elements, dtype=str)
        self.atoms    = np.ones_like(self.elements,        dtype=int)
        #
        ### Harmonic Frequencies
        #
        if harmonic: frequencies = harmonic["freq_wavenumber"]
        #
        ### Rotational data
        #
        self.rotationalsymm = 1.0
        self.rotationaltemp = np.pi/np.prod( self._rotation_constant() )

        if beta_orbital_occup.size:
            fermi_energy = max(alpha_orbital_energy[alpha_orbital_occup==1.0][-1], beta_orbital_energy[beta_orbital_occup==1.0][-1])
            self.orbitalenergies = np.concatenate((alpha_orbital_energy, beta_orbital_energy)) - fermi_energy

        else:
            fermi_energy = alpha_orbital_energy[alpha_orbital_occup==1.0][-1]
            self.orbitalenergies = np.concatenate((alpha_orbital_energy, alpha_orbital_energy)) - fermi_energy

        if not frequencies: frequencies = 3*[0.0]
        if any(i < 0 for i in frequencies): raise RuntimeError(f"found imaginary frequency")

        frequencies = np.asarray(frequencies, dtype=float)
        frequencies = frequencies[frequencies != 0.0]
        
        self.frequencies     = frequencies*cm1_to_hz
        self.zeropointenergy = 0.5*np.sum(frequencies)*cm1_to_ev

        self.volume = 0.0
    
        return multiplicity
#
### Set magnetic moments per atom
#
    def _set_magnetic_moments(self, multiplicity):

        if self.magnetization is not None: return self.magnetization

        if multiplicity <= 0.0:
            raise RuntimeError(f"cannot resolve magnetic moments. Use 'magnetization=' instead")

        magnetization = np.array([multiplicity-1.0 if i in spin_crossover_metals else 0.0 for i in self.elements], dtype=float)

        return magnetization
#
### Determine rotational constant
#
    def _rotation_constant(self):

        mass = self.meanfield.mol.atom_mass_list(isotope_avg=True)
        coordinates = self.meanfield.mol.atom_coords()

        center_of_mass = np.einsum("z,zr->r", mass, coordinates) / mass.sum()

        distance = coordinates - center_of_mass
        moment_of_inertia = np.einsum("z,zr,zs->rs", mass, distance, distance)
        moment_of_inertia = np.eye(3)*moment_of_inertia.trace() - moment_of_inertia

        rotation_constant = np.sort( np.linalg.eigvalsh(moment_of_inertia) )

        unit_hz = h_bar/(4.0*np.pi*atomic_mass*bohr_radius*bohr_radius)

        with np.errstate(divide="ignore"):
            result = h_over_k*unit_hz/rotation_constant
            
        return result

#############################################################################################################################

##     ##  #######  ########  ##     ## ##       ########  ######
###   ### ##     ## ##     ## ##     ## ##       ##       ##    ##
#### #### ##     ## ##     ## ##     ## ##       ##       ##
## ### ## ##     ## ##     ## ##     ## ##       ######    ######
##     ## ##     ## ##     ## ##     ## ##       ##             ##
##     ## ##     ## ##     ## ##     ## ##       ##       ##    ##
##     ##  #######  ########   #######  ######## ########  ######

#############################################################################################################################
### IO operations

def fileio(filename=None, option="read"):

    if option in ("w", "write"):

        try:
            return True if os.path.getsize(filename) == 0 else False

        except OSError:
            return True

    elif option in ("r","s","read","silent"):

        try:
            if os.path.getsize(filename) == 0:
                if option in ("r","read"):
                    raise RuntimeError(f"{filename} is empty")
                return False
            else:
                return True

        except OSError:
            if option in ("r","read"):
                raise FileNotFoundError(f"{filename} not found")
            return False
    
#############################################################################################################################
