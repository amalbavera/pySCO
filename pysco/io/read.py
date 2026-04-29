#############################################################################################################################

##       #### ########  ########     ###    ########  #### ########  ######
##        ##  ##     ## ##     ##   ## ##   ##     ##  ##  ##       ##    ##
##        ##  ##     ## ##     ##  ##   ##  ##     ##  ##  ##       ##
##        ##  ########  ########  ##     ## ########   ##  ######    ######
##        ##  ##     ## ##   ##   ######### ##   ##    ##  ##             ##
##        ##  ##     ## ##    ##  ##     ## ##    ##   ##  ##       ##    ##
######## #### ########  ##     ## ##     ## ##     ## #### ########  ######

#############################################################################################################################

import numpy as np

from pathlib import Path as path
#
### local libraries
#
from ..utils.constants import spin_crossover_metals
from ..utils.constants import cm1_to_hz, cm1_to_ev, au_to_ev
from ..utils.constants import h_over_k, c_light, h_bar, bohr_radius, atomic_mass

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

from .files   import fileio

#############################################################################################################################
### ATTENTION: Information regarding units from output files

# spin              = > up-electrons - down-electrons
# orbit             = > dimensionless
# atoms             = > dimensionless 
# energy            = > eV
# volume            = > Angstroem^3
# elements          = > dimensionless
# jahnteller        = > dimensionless
# frequencies       = > Hz
# configuration     = > dimensionless
# orbital_energies  = > eV
# zero_point_energy = > eV

# rotational_symmetry    = > dimensionless
# rotational_temperature => K

#############################################################################################################################

########     ###     ######  ######## ##       #### ##    ## ########
##     ##   ## ##   ##    ## ##       ##        ##  ###   ## ##
##     ##  ##   ##  ##       ##       ##        ##  ####  ## ##
########  ##     ##  ######  ######   ##        ##  ## ## ## ######
##     ## #########       ## ##       ##        ##  ##  #### ##
##     ## ##     ## ##    ## ##       ##        ##  ##   ### ##
########  ##     ##  ######  ######## ######## #### ##    ## ########

#############################################################################################################################

class reader:
    """
    Provides the base-line class with the purpose of reading the outputs from electronic structure calculations.

    Parameters
    ----------
        jahnteller : int, default=1
            Number of Jahn-Teller distortions for the given spin state.
        
        orbit : int, default=0
            Angular momentum number, L, for use in the calculation of the orbital contribution to the spin entropy.

        configs : int, default=1
            Number of symmetry related spin configurations possible in intermediate spin states for the calculation of the
            configuration entropy.

    Modules
    -------
        set_spin
            Set the spin moment per atom from a user-defined list.

        set_orbit
            Set the orbital angular momentum per atom from a user-defined list.

    Private Modules
    ---------------
        _set_magnetic_moments
            Prepare the list with the spin moment per atom.

        _set_orbital_momentum
            Prepare the list with the orbital angular momentum per atom.

        _set_zero_point_energy
            Compute the zero-point energy from the extracted harmonic frequencies.

        _check_attributes
            Verify that all attributes were properly set.

    Attributes
    ----------
        atoms : ndarray[tuple[str], dtype[str]]
            Total count for each of the unique elements in the system.

        energy : float, eV
            Self-consistent total energy.

        volume : float, Angstroem^3
            Volume of the unit cell.

        elements : ndarray[tuple[int], dtype[int]]
            Unique elements in the system.

        frequencies : ndarray[tuple[float], dtype[float]], Hz
            Harmonic vibrational frequencies.

        orbital_energies : ndarray[tuple[float], dtype[float]], eV
            Single-particle energies subtracted by the Fermi energy.

        zero_point_energy : float, eV
            Zero-point energy computed from the set of harmonic vibrational frequencies.

        rotational_symmetry : float, dimensionless
            Total number of symmetry-related rotation operations for the system.

        rotational_temperature : float, K
            Rotational temperature obtained from the moments of inertia of the system.
    """

    __slots__ = (
                 "spin", "atoms", "orbit", "energy", "volume", "elements", "rotational_symmetry", "rotational_temperature",
                 "jahnteller", "frequencies", "configurations", "orbital_energies", "zero_point_energy"
                )
#
### Initialize
#
    def __init__(self, jahnteller:int=1, orbit:int=0, configs:int=1) -> None:

        self.spin              = None
        self.orbit             = None
        self.atoms             = None
        self.energy            = None
        self.volume            = None
        self.elements          = None
        self.jahnteller        = None
        self.frequencies       = None
        self.configurations    = None
        self.orbital_energies  = None
        self.zero_point_energy = None

        self.rotational_symmetry    = 0.0
        self.rotational_temperature = 0.0

        try:
            orbit = float(orbit)
            if orbit < 0:
                raise RuntimeError("orbit must be a positive number.")
        except:
            raise RuntimeError(f"orbit = {orbit} is not a number.")

        try:
            self.jahnteller = float(jahnteller)
            if jahnteller < 1:
                raise RuntimeError("jahnteller must be a positive non-zero number.")
        except:
            raise RuntimeError(f"jahnteller = {jahnteller} is not a number.")

        try:
            self.configurations = float(configs)
            if configs < 1:
                raise RuntimeError("configs must be a positive non-zero number.")
        except:
            raise RuntimeError(f"configs = {configs} is not a number.")
#
### PUBLIC: Define magnetic moments per atom
#
    def set_spin(self, spin:list[int]) -> None:
        """
        Set the spin moment per atom from a user-defined list.

        Parameters
        ----------
            spin : list[int]
                List containing the spin moment per atom.

        Returns
        -------
            None
        """

        if isinstance(spin, list): spin = np.asarray(spin, dtype=float)

        if not spin.size == self.atoms.size:
            raise RuntimeError(f"size missmatch between spin ({spin.size}) and atoms ({self.atoms.size})")

        self.spin = spin
#
### PUBLIC: Define orbital momentum per atom
#
    def set_orbit(self, orbit:list[int]) -> None:
        """
        Set the orbital angular momentum per atom from a user-defined list.

        Parameters
        ----------
            orbit : list[int]
                List containing the orbital angular momentum per atom.

        Returns
        -------
            None
        """

        if isinstance(orbit, list): orbit = np.asarray(orbit, dtype=float)

        if not orbit.size == self.atoms.size:
            raise RuntimeError(f"size missmatch between orbit ({orbit.size}) and atoms ({self.atoms.size})")

        self.orbit = orbit
#
### PRIVATE: Set magnetic moments per atom
#
    def _set_magnetic_moments(self, multiplicity:float) -> np.ndarray[tuple[float], np.dtype[float]]:
        """
        Prepare the list with the spin moment per atom. The spin moment is assigned by comparing each element
        with respect to the predefiend tuple ('Cr', 'Mn', 'Fe', 'Co'). A defaul value of zero is
        used for the rest of the elements

        Parameters
        ----------
            multiplicity : float
                Multiplicity for the spin-crossover metal center.

        Returns
        -------
            spin : ndarray[tuple[float], dtype[float]]
        """

        if multiplicity == 0.0:
            raise RuntimeError(f"cannot resolve magnetic moments. Use 'spin=' instead.")
        
        atoms = np.repeat(self.elements, self.atoms)

        spin  = np.array([multiplicity-1.0 if i in spin_crossover_metals else 0.0 for i in atoms], dtype=float)

        return spin
#
### PRIVATE: Set orbital momentum per atom
#
    def _set_orbital_momentum(self, momentum:float) -> np.ndarray[tuple[float], np.dtype[float]]:
        """
        Prepare the list with the orbital angular momentum per atom. The orbital angular momentum is assigned by
        comparing each element with respect to the predefiend tuple ('Cr', 'Mn', 'Fe', 'Co'). A defaul value of
        zero is used for the rest of the elements

        Parameters
        ----------
            momentum : float
                Angular orbital momentum for the spin-crossover metal center.

        Returns
        -------
            orbit : ndarray[tuple[float], dtype[float]]
        """

        atoms = np.repeat(self.elements, self.atoms)

        orbit = np.array([momentum if i in spin_crossover_metals else 0.0 for i in atoms], dtype=float)

        return orbit
#
### PRIVATE: Compute the zero-point energy from harmonic frequencies
#
    def _set_zero_point_energy(self) -> float:
        """
        Compute the zero-point energy from the extracted harmonic frequencies.

        Parameters
        ----------
            None

        Returns
        -------
            None
        """

        zero_point_energy = 0.5*np.sum(self.frequencies)*cm1_to_ev/cm1_to_hz

        return zero_point_energy
#
### PRIVATE: Sanity check for all mandatory attributes
#
    def _check_attributes(self) -> None:
        """
        Sanity check to verify no unassigned attributes remaining that could cause errors during runtime.

        Parameters
        ----------
            None

        Returns
        -------
            None
        """

        for attribute in self.__slots__:

            if self.__getattribute__(attribute) is None:
                raise RuntimeError(f"attribute {attribute} must not be None.")
            
            if isinstance(self.__getattribute__(attribute), list):
                self.__setattr__(attribute, np.asarray(self.__getattribute__(attribute)))

#############################################################################################################################

##     ##    ###     ######  ########  
##     ##   ## ##   ##    ## ##     ## 
##     ##  ##   ##  ##       ##     ## 
##     ## ##     ##  ######  ########  
 ##   ##  #########       ## ##        
  ## ##   ##     ## ##    ## ##        
   ###    ##     ##  ######  ##

#############################################################################################################################

class Vasp(reader):
    """
    Read outputs from the Vienna Ab initio Simulation Package.

    Parameters
    ----------
        filesdir : str|path
            Location of the directory containing the OUTCAR, POSCAR, EIGENVAL, MAGCAR, and VIBCAR

        jahnteller : int, default=1
            Number of Jahn-Teller distortions for the given spin state.
        
        orbit : int, default=0
            Angular momentum number, L, for use in the calculation of the orbital contribution to the spin entropy.

        configs : int, default=1
            Number of symmetry related spin configurations possible in intermediate spin states for the calculation of the
            configuration entropy.

    Modules
    -------
        set_spin
            Set the spin moment per atom from a user-defined list.

        set_orbit
            Set the orbital angular momentum per atom from a user-defined list.

    MAGCAR
    ------
        Contains the analogous magnetization as reported in the OUTCAR, but with the correct rounded spin moments, e.g.,
        ```
        magnetization (x)

        # of ion       s       p       d       tot
        ------------------------------------------
            1        0.000   0.000   4.000   4.000
            2        0.000   0.000   4.000   4.000
            3        0.000   0.000   0.000   0.000
            ...
        256        0.000   0.000   0.000   0.000
        ------------------------------------------
        tot          0.000   0.000   8.000   8.000
        ```

    VIBCAR
    ------
        Is the OUTCAR that contains the harmonic vibrational frequencies. Either rename the original OUTCAR or save time by
        running the following command from the console,
        ```
        grep 2PiTHz OUTCAR > VIBCAR
        ```
    """
#
### Initialize
#
    def __init__(self, filesdir:str|path, jahnteller:int=1, orbit:int=0, configs:int=1) -> None:
        super().__init__(jahnteller=jahnteller, orbit=orbit, configs=configs)

        if isinstance(filesdir,str): filesdir = path(filesdir)

        spin       = self._get_output_data(filesdir)

        self.spin  = self._set_magnetic_moments(spin)
        self.orbit = self._set_orbital_momentum(orbit)

        self.zero_point_energy = self._set_zero_point_energy()

        self._check_attributes()
#
### PRIVATE: Extract all data from output directory
#
    def _get_output_data(self, filesdir:path) -> list:

        outfile = filesdir/"OUTCAR"
        posfile = filesdir/"POSCAR"
        eigfile = filesdir/"EIGENVAL"

        magfile = filesdir/"MAGCAR"
        vibfile = filesdir/"VIBCAR"

        fileio(outfile, option="r")
        fileio(posfile, option="r")
        fileio(eigfile, option="r")

        skipped_lines = 0

        spin        = []
        frequencies = []
        #
        ### Runtime variables
        #
        found_spin     = False
        skip_mag_lines = 3
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
        nspin    = int(data[0].rstrip().split()[-1])
        kpoints  = int(details[1])
        bands    = int(details[2])

        points   = [j for i in range(8,len(data), bands+2) for j in range(i,i+bands)]

        alpha    = np.loadtxt(data[points], usecols=1, dtype=float)

        if nspin == 2: beta = np.loadtxt(data[points], usecols=2, dtype=float)

        self.orbital_energies = alpha - fermi_energy if nspin == 1 else np.concatenate((alpha,beta)) - fermi_energy
        #
        ### Harmonic Frequencies
        #
        if fileio(vibfile, option="s"):
            with open(vibfile, "r") as out:
                for line in out:
                    line = line.strip()

                    if vasp_frequency_str in line:
                        if vasp_imag_freq_str in line: raise RuntimeError(f"found imaginary frequency in {vibfile.name}.")
                        frequencies += [float(line.split()[7])]
        #
        ### Spin
        #
        if (self.spin is None) and fileio(magfile, option="s"):
            with open(magfile, "r") as out:
                for line in out:
                    line = line.strip()

                    if vasp_mag_start_str in line:
                        skipped_lines = 0
                        found_spin = True
                        continue

                    if   (found_spin) and (skipped_lines < skip_mag_lines):
                        skipped_lines += 1
                        continue
                    elif (found_spin) and (vasp_mag_stop_str in line):
                        skipped_lines = 0
                        found_spin = False
                        continue
                    elif (found_spin) and (skipped_lines == skip_mag_lines):
                        spin += [ float(line.split()[-1]) ]
                        continue

        if not frequencies: frequencies = 3*[0.0]

        frequencies = np.asarray(frequencies, dtype=float)
        frequencies = frequencies[frequencies != 0.0]

        self.frequencies = frequencies*cm1_to_hz

        return spin
#
### PRIVATE: Set magnetic moments
#
    def _set_magnetic_moments(self, spin:list) -> list:

        if not spin: raise RuntimeError("please provide a MAGCAR.")

        spin = np.asarray(spin, dtype=float)

        return spin
        
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

class Gaussian(reader):
    """
    Read output from Gaussian.

    Parameters
    ----------
        outfile : str|path
            Location of the output file.

        jahnteller : int, default=1
            Number of Jahn-Teller distortions for the given spin state.
        
        orbit : int, default=0
            Angular momentum number, L, for use in the calculation of the orbital contribution to the spin entropy.

        configs : int, default=1
            Number of symmetry related spin configurations possible in intermediate spin states for the calculation of the
            configuration entropy.

    Modules
    -------
        set_spin
            Set the spin moment per atom from a user-defined list.

        set_orbit
            Set the orbital angular momentum per atom from a user-defined list.
    """
#
### Initialize
#
    def __init__(self, outfile:str|path, jahnteller:int=1, orbit:int=0, configs:int=1) -> None:
        super().__init__(jahnteller=jahnteller, orbit=orbit, configs=configs)

        if isinstance(outfile,str): outfile = path(outfile)
        fileio(outfile, option="r")

        multiplicity = self._get_output_data(outfile)

        self.spin    = self._set_magnetic_moments(multiplicity)
        self.orbit   = self._set_orbital_momentum(orbit)

        self.zero_point_energy = self._set_zero_point_energy()

        self._check_attributes()
#
### PRIVATE: Extract all data from output file
#
    def _get_output_data(self, data:path) -> float:

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
                    self.rotational_symmetry = 1.0/float(line.split()[3])

                if gaussian_rot_temp_str in line:
                    self.rotational_temperature = np.pi/np.prod([ float(i) for i in line.split()[3:] ])

        alpha_orbital_energy = np.asarray(alpha_orbital_energy, dtype=float)*au_to_ev
        beta_orbital_energy  = np.asarray(beta_orbital_energy,  dtype=float)*au_to_ev
        alpha_orbital_occup  = np.asarray(alpha_orbital_occup,  dtype=float)
        beta_orbital_occup   = np.asarray(beta_orbital_occup,   dtype=float)

        if beta_orbital_occup.size:
            fermi_energy = max(alpha_orbital_energy[alpha_orbital_occup==1.0][-1], beta_orbital_energy[beta_orbital_occup==1.0][-1])
            self.orbital_energies = np.concatenate((alpha_orbital_energy, beta_orbital_energy)) - fermi_energy

        else:
            fermi_energy = alpha_orbital_energy[alpha_orbital_occup==1.0][-1]
            self.orbital_energies = np.concatenate((alpha_orbital_energy, alpha_orbital_energy)) - fermi_energy

        self.elements = np.asarray(elements,   dtype=str)
        self.atoms    = np.ones_like(elements, dtype=int)

        if not frequencies: frequencies = 3*[0.0]
        if any(i < 0 for i in frequencies): raise RuntimeError(f"found imaginary frequency in {data.name}.")

        frequencies = np.asarray(frequencies, dtype=float)
        frequencies = frequencies[frequencies != 0.0]

        self.frequencies = frequencies*cm1_to_hz

        self.volume = 0.0

        return multiplicity

#############################################################################################################################

##    ## ##      ##  ######  ##     ## ######## ##     ## 
###   ## ##  ##  ## ##    ## ##     ## ##       ###   ### 
####  ## ##  ##  ## ##       ##     ## ##       #### #### 
## ## ## ##  ##  ## ##       ######### ######   ## ### ## 
##  #### ##  ##  ## ##       ##     ## ##       ##     ## 
##   ### ##  ##  ## ##    ## ##     ## ##       ##     ## 
##    ##  ###  ###   ######  ##     ## ######## ##     ##

#############################################################################################################################

class NWChem(reader):
    """
    Read output from NWChem.

    Parameters
    ----------
        outfile : str|path
            Location of the output file.

        jahnteller : int, default=1
            Number of Jahn-Teller distortions for the given spin state.
        
        orbit : int, default=0
            Angular momentum number, L, for use in the calculation of the orbital contribution to the spin entropy.

        configs : int, default=1
            Number of symmetry related spin configurations possible in intermediate spin states for the calculation of the
            configuration entropy.

    Modules
    -------
        set_spin
            Set the spin moment per atom from a user-defined list.

        set_orbit
            Set the orbital angular momentum per atom from a user-defined list.
    """
#
### Initialize
#
    def __init__(self, outfile:str|path, jahnteller:int=1, orbit:int=0, configs:int=1) -> None:
        super().__init__(jahnteller=jahnteller, orbit=orbit, configs=configs)

        if isinstance(outfile,str): outfile = path(outfile)
        fileio(outfile, option="r")

        multiplicity = self._get_output_data(outfile)

        self.spin    = self._set_magnetic_moments(multiplicity)
        self.orbit   = self._set_orbital_momentum(orbit)

        self.zero_point_energy = self._set_zero_point_energy()

        self._check_attributes()
#
### PRIVATE: Extract all data from output file
#
    def _get_output_data(self, data:path) -> float:

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
                    self.rotational_symmetry = 1.0/float(line.split()[8].replace(")",""))
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
                    self.rotational_temperature = np.pi/np.prod(rotational_temperatures)
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
            self.orbital_energies = np.concatenate((alpha_orbital_energy, beta_orbital_energy)) - fermi_energy

        else:
            fermi_energy = alpha_orbital_energy[alpha_orbital_occup==1.0][-1]
            self.orbital_energies = np.concatenate((alpha_orbital_energy, alpha_orbital_energy)) - fermi_energy

        self.elements = np.asarray(elements,   dtype=str)
        self.atoms    = np.ones_like(elements, dtype=int)

        if not frequencies: frequencies = 3*[0.0]
        if any(i < 0 for i in frequencies): raise RuntimeError(f"found imaginary frequency in {data.name}.")

        frequencies = np.asarray(frequencies, dtype=float)
        frequencies = frequencies[frequencies != 0.0]

        self.frequencies = frequencies*cm1_to_hz

        self.volume = 0.0

        return multiplicity

#############################################################################################################################

 #######  ########   ######     ###    
##     ## ##     ## ##    ##   ## ##   
##     ## ##     ## ##        ##   ##  
##     ## ########  ##       ##     ## 
##     ## ##   ##   ##       ######### 
##     ## ##    ##  ##    ## ##     ## 
 #######  ##     ##  ######  ##     ##

#############################################################################################################################

class Orca(reader):
    """
    Read output from Orca.

    Parameters
    ----------
        outfile : str|path
            Location of the output file.

        jahnteller : int, default=1
            Number of Jahn-Teller distortions for the given spin state.
        
        orbit : int, default=0
            Angular momentum number, L, for use in the calculation of the orbital contribution to the spin entropy.

        configs : int, default=1
            Number of symmetry related spin configurations possible in intermediate spin states for the calculation of the
            configuration entropy.

    Modules
    -------
        set_spin
            Set the spin moment per atom from a user-defined list.

        set_orbit
            Set the orbital angular momentum per atom from a user-defined list.
    """
#
### Initialize
#
    def __init__(self, outfile:str|path, jahnteller:int=1, orbit:int=0, configs:int=1) -> None:
        super().__init__(jahnteller=jahnteller, orbit=orbit, configs=configs)

        if isinstance(outfile,str): outfile = path(outfile)
        fileio(outfile, option="r")

        multiplicity = self._get_output_data(outfile)

        self.spin    = self._set_magnetic_moments(multiplicity)
        self.orbit   = self._set_orbital_momentum(orbit)

        self.zero_point_energy = self._set_zero_point_energy()

        self._check_attributes()
#
### PRIVATE: Extract all data from output file
#
    def _get_output_data(self, data:path) -> float:

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
                    self.rotational_symmetry = 1.0/float(line.split()[5])

                if orca_rot_const_str in line:
                    rotational_temperatures     = np.asarray(line.split()[4:], dtype=float)*1e+02*h_over_k*c_light
                    self.rotational_temperature = np.pi/np.prod(rotational_temperatures)

        alpha_orbital_energy = np.asarray(alpha_orbital_energy, dtype=float)
        beta_orbital_energy  = np.asarray(beta_orbital_energy,  dtype=float)
        alpha_orbital_occup  = np.asarray(alpha_orbital_occup,  dtype=float)
        beta_orbital_occup   = np.asarray(beta_orbital_occup,   dtype=float)

        if beta_orbital_occup.size:
            fermi_energy = max(alpha_orbital_energy[alpha_orbital_occup==1.0][-1], beta_orbital_energy[beta_orbital_occup==1.0][-1])
            self.orbital_energies = np.concatenate((alpha_orbital_energy, beta_orbital_energy)) - fermi_energy

        else:
            fermi_energy = alpha_orbital_energy[alpha_orbital_occup==1.0][-1]
            self.orbital_energies = np.concatenate((alpha_orbital_energy, alpha_orbital_energy)) - fermi_energy

        self.elements = np.asarray(elements,   dtype=str)
        self.atoms    = np.ones_like(elements, dtype=int)

        if not frequencies: frequencies = 3*[0.0]
        if any(i < 0 for i in frequencies): raise RuntimeError(f"found imaginary frequency in {data.name}.")

        frequencies = np.asarray(frequencies, dtype=float)
        frequencies = frequencies[frequencies != 0.0]

        self.frequencies = frequencies*cm1_to_hz

        self.volume = 0.0

        return multiplicity

#############################################################################################################################

########  ##    ##  ######   ######  ########
##     ##  ##  ##  ##    ## ##    ## ##
##     ##   ####   ##       ##       ##
########     ##     ######  ##       ######
##           ##          ## ##       ##
##           ##    ##    ## ##    ## ##
##           ##     ######   ######  ##

#############################################################################################################################

class pySCF(reader):
    """
    Read output from pySCF.

    Parameters
    ----------
        meanfield : SCF object
            Mean-field object generated with pySCF.

        harmonic_analysis : dict
            Dictionary created through the `thermo` module in pySCF.

        jahnteller : int, default=1
            Number of Jahn-Teller distortions for the given spin state.
        
        orbit : int, default=0
            Angular momentum number, L, for use in the calculation of the orbital contribution to the spin entropy.

        configs : int, default=1
            Number of symmetry related spin configurations possible in intermediate spin states for the calculation of the
            configuration entropy.

    Modules
    -------
        set_spin
            Set the spin moment per atom from a user-defined list.

        set_orbit
            Set the orbital angular momentum per atom from a user-defined list.

    Extra
    -----
        meanfield
            Object created with pySCF. It often is created from the `pyscf.dft` module, for instance,
        ```
        from pyscf import gto, dft

        mol = gto.M('file.xyz')
        meanfield = dft.UKS(mol)
        ```
        harmonic_analysis
            Output dictionary from the `thermo` module in pySCF. It results from running an analogous code to
        ```
        from pyscf.hessian import thermo

        hessian = meanfield.Hessian().kernel()
        harmonic_analysis = thermo.thermo(meanfield, hessian)
        ```
    """
#
### Initialize
#
    def __init__(self, meanfield, harmonic_analysis:dict=None, jahnteller:int=1, orbit:int=0, configs:int=1) -> None:
        super().__init__(jahnteller=jahnteller, orbit=orbit, configs=configs)

        multiplicity = self._get_output_data(meanfield, harmonic_analysis)

        self.spin    = self._set_magnetic_moments(multiplicity)
        self.orbit   = self._set_orbital_momentum(orbit)

        self.zero_point_energy = self._set_zero_point_energy()

        self._check_attributes()
#
### PRIVATE: Extract all data from output objects
#
    def _get_output_data(self, meanfield, harmonic) -> float:

        alpha_orbital_energy = []
        beta_orbital_energy  = []
        alpha_orbital_occup  = []
        beta_orbital_occup   = []

        frequencies          = []
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
        self.rotational_symmetry = 1.0
        self.rotational_temperature = np.pi/np.prod( self._rotation_constant(meanfield) )

        if beta_orbital_occup.size:
            fermi_energy = max(alpha_orbital_energy[alpha_orbital_occup==1.0][-1], beta_orbital_energy[beta_orbital_occup==1.0][-1])
            self.orbital_energies = np.concatenate((alpha_orbital_energy, beta_orbital_energy)) - fermi_energy

        else:
            fermi_energy = alpha_orbital_energy[alpha_orbital_occup==1.0][-1]
            self.orbital_energies = np.concatenate((alpha_orbital_energy, alpha_orbital_energy)) - fermi_energy

        if not frequencies: frequencies = 3*[0.0]
        if any(i < 0 for i in frequencies): raise RuntimeError(f"found imaginary frequency.")

        frequencies = np.asarray(frequencies, dtype=float)
        frequencies = frequencies[frequencies != 0.0]

        self.frequencies = frequencies*cm1_to_hz

        self.volume = 0.0

        return multiplicity
#
### Determine rotational constant
#
    @staticmethod
    def _rotation_constant(meanfield) -> float:

        mass = meanfield.mol.atom_mass_list(isotope_avg=True)
        coordinates = meanfield.mol.atom_coords()

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
