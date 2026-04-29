#############################################################################################################################

##       #### ########  ########     ###    ########  #### ########  ######
##        ##  ##     ## ##     ##   ## ##   ##     ##  ##  ##       ##    ##
##        ##  ##     ## ##     ##  ##   ##  ##     ##  ##  ##       ##
##        ##  ########  ########  ##     ## ########   ##  ######    ######
##        ##  ##     ## ##   ##   ######### ##   ##    ##  ##             ##
##        ##  ##     ## ##    ##  ##     ## ##    ##   ##  ##       ##    ##
######## #### ########  ##     ## ##     ## ##     ## #### ########  ######

#############################################################################################################################

import math
import numpy as np
#
### local libraries
#
from ..io.read         import reader
from ..utils.constants import spin_crossover_metals
from ..utils.constants import ev_to_joule, ang3_to_m3
from ..utils.constants import k_times_n, atm_pressure, h_over_k
from ..utils.constants import exp_thresh, log_thresh

#############################################################################################################################

######## ##     ## ######## ########  ##     ##  #######  ########  ##    ## ##    ##    ###    ##     ## ####  ######  
   ##    ##     ## ##       ##     ## ###   ### ##     ## ##     ##  ##  ##  ###   ##   ## ##   ###   ###  ##  ##    ## 
   ##    ##     ## ##       ##     ## #### #### ##     ## ##     ##   ####   ####  ##  ##   ##  #### ####  ##  ##       
   ##    ######### ######   ########  ## ### ## ##     ## ##     ##    ##    ## ## ## ##     ## ## ### ##  ##  ##       
   ##    ##     ## ##       ##   ##   ##     ## ##     ## ##     ##    ##    ##  #### ######### ##     ##  ##  ##       
   ##    ##     ## ##       ##    ##  ##     ## ##     ## ##     ##    ##    ##   ### ##     ## ##     ##  ##  ##    ## 
   ##    ##     ## ######## ##     ## ##     ##  #######  ########     ##    ##    ## ##     ## ##     ## ####  ######

#############################################################################################################################

class Thermodynamic:
    """
    Provides the base-line class for a thermodynamic spin-crossover model.

    Parameters
    ----------
        ls : reader
            Defines the reference low-spin state read from an electronic structure output.

        hs : reader
            Defines the reference high-spin state read from an electronic structure output.

        centers : float, default=1.0
            Number of spin-crossover metal centers used to normalize reported data.

        metals : tuple(str), default=('Cr', 'Mn', 'Fe', 'Co')
            Define the transition metals of interest for automatic determination of the parameter `centers`.

    Modules
    -------
        set_nhs_list
            Define a list of points for the relative high-spin population, or magnetization, within the interval
            `min_nhs` <= n_HS <= `max_nhs`.
        
        set_min_nhs
            Set the lower threshold for the relative high-spin population.

        set_max_nhs
            Set the upper threshold for the relative high-spin population.

    Attributes
    ----------
        min_nhs : float, default=1.00E-06
            Lower bound for the relative high-spin population of interest.

        max_nhs : float, default=9.99E-01
            Upper bound for the relative high-spin population of interest.
    """

    __slots__ = (
                 "ls", "hs", "atoms", "centers",
                 "PdV",
                 "min_nhs", "max_nhs"
                )
#
### Initialize
#
    def __init__(self, ls:reader, hs:reader, centers:float=1.0, metals:tuple[str]=spin_crossover_metals) -> None:

        diff_str = "in low- and high-spin differ in size."

        assert all(ls.elements == hs.elements), f"elements {diff_str}"
        assert all(ls.atoms    == hs.atoms),    f"atoms {diff_str}"

        assert ls.frequencies.shape      == hs.frequencies.shape,      f"harmonic frequencies {diff_str}"
        assert ls.orbital_energies.shape == hs.orbital_energies.shape, f"single-particle energies {diff_str}"

        self.ls = ls
        self.hs = hs

        self.atoms   = np.repeat(ls.elements, ls.atoms)

        self.centers = np.count_nonzero( np.isin(self.atoms[ls.spin != hs.spin], metals) ) if not centers else centers

        self.PdV = self._PdV()

        self.min_nhs = 1.00e-06
        self.max_nhs = 9.99e-01
#
### PUBLIC: List of relative high-spin population points
#
    @staticmethod
    def set_nhs_list(points:int, min_nhs:float=1.00e-06, max_nhs:float=9.99e-01) -> np.ndarray[tuple[float], np.dtype[float]]:
        """
        Define a list of points for the relative high-spin population, or magnetization, n_HS within the interval
        `min_nhs` <= n_HS <= `max_nhs`.

        Parameters
        ----------
            points : int
                Total number of point for the list.

            min_nhs : float, default=1.00E-06
                Lower bound for for the relative high-spin population.

            max_nhs : float, default=9.99E-01
                Upper bound for the relative high-spin population.

        Returns
        -------
            1D array : ndarray[tuple[float], dtype[float]]

        Units
        -----
            Dimensionless
        """

        lower_bound = -math.log(1.0/min_nhs - min_nhs)
        upper_bound = -math.log(1.0/max_nhs - max_nhs)
        
        points_list = np.linspace(lower_bound, upper_bound, num=points)

        nhs_points  = 1.0/(1.0 + np.exp(-points_list))

        return nhs_points
#
### PUBLIC: Set lower limit for the relative high-spin population
#
    def set_min_nhs(self, nhs:float) -> None:
        """
        Set the lower limit allowed for the relative high-spin population.

        Parameters
        ----------
            nhs : float
                Define the lower limit subject to 0 < nhs < `max_nhs`.

        Returns
        -------
            None
        """

        if nhs == 0.0:
            raise RuntimeWarning("a value of 0.0 will likely cause numerical instabilities.")
        
        if nhs <  0.0:
            raise RuntimeError("the relative high-spin population cannot be negative.")
        
        if nhs >= self.max_nhs:
            raise RuntimeError("the lower limit cannot be larger than the upper limit.")
        
        if not (0.0 < nhs < self.max_nhs):
            raise RuntimeError("lower limit out of bounds.")

        self.min_nhs = nhs

### PUBLIC: Set lower limit for the relative high-spin population
#
    def set_max_nhs(self, nhs:float) -> None:
        """
        Set the upper limit allowed for the relative high-spin population.

        Parameters
        ----------
            nhs : float
                Define the upper limit subject to `min_nhs` < nhs < 1.

        Returns
        -------
            None
        """

        if nhs == 1.0:
            raise RuntimeWarning("a value of 1.0 will likely cause numerical instabilities.")
        
        if nhs <  0.0:
            raise RuntimeError("the relative high-spin population cannot be negative.")
        
        if nhs <= self.min_nhs:
            raise RuntimeError("the upper limit cannot be smaller than the lower limit.")
        
        if not (self.min_nhs < nhs < 1.0):
            raise RuntimeError("upper limit out of bounds.")

        self.max_nhs = nhs
#
### PUBLIC: Spin-crossover energy gap
#
    def spin_gap(self, spin_state:reader=None, zero_point_energy:bool=False, convert:float=ev_to_joule) -> float:
        """
        Compute the spin-crossover energy, between a high- and the low-spin state normalized to the number of
        spin transition centers.

        Parameters
        ----------
            spin_state : reader. default=high-spin
                High-spin state of interest. May be an intermediate spin state.

            zero_point_energy : bool, default=False
                Add the zero-point correction to the reported value.

            convert : float, default=ev_to_joule
                Conversion factor from eV to desired units.

        Returns
        -------
            dE_ele : float
                Spin-crossover energy gap

        Units
        -----
            J/mol
        """

        low_spin  = self.ls
        high_spin = self.hs if spin_state is None else spin_state

        dE_ele = (high_spin.energy - low_spin.energy)*convert

        if zero_point_energy:
            dE_ele += (high_spin.zero_point_energy - low_spin.zero_point_energy)*convert

        return dE_ele
#
### PULBIC: Compute change in Gibbs free energy
#
    def delta_gibbs(self, T:np.ndarray[tuple[float], np.dtype[float]]) -> tuple[np.ndarray[tuple[float], np.dtype[float]], np.ndarray[tuple[float], np.dtype[float]], np.ndarray[tuple[float], np.dtype[float]]]:
        """
        Compute the change in the Gibbs free energy between the low- and high-spin states.

        ΔG = G_HS - G_LS = dH - T\*ΔS, where

        ΔH = ΔE_ele + ΔE_vib + ΔE_rot + P\*ΔV

        ΔS = ΔS_ele + ΔS_vib + ΔS_rot + ΔS_config

        with ΔS_ele = ΔS_spin + ΔS_orb + ΔS_Fermi

        Parameters
        ----------
            T : ndarray[tuple[float], dtype[float]]
                Temperature(s) of interest in units of K.

        Returns
        ------
            dG : ndarray[tuple[float], dtype[float]]
                Gibbs free energy change.

            dH : ndarray[tuple[float], dtype[float]]
                Enthalpy change.

            dS : ndarray[tuple[float], dtype[float]]
                Entropy change.

        Units
        -----
            J/mol, J/mol, J/(mol K), respectivel
        """

        if (T.size == 1) and ( (T <= 0.0) or (T is np.nan) ):
            inf = np.array([-float("inf")])

            return inf, inf, inf

        E_vib_ls, S_vib_ls = self._E_vib_and_S_vib(T, self.ls)
        E_vib_hs, S_vib_hs = self._E_vib_and_S_vib(T, self.hs)

        E_ele_ls  = self.ls.energy*ev_to_joule
        E_ele_hs  = self.hs.energy*ev_to_joule

        S_ele_ls  = self._S_ele(T, self.ls)
        S_ele_hs  = self._S_ele(T, self.hs)

        S_rot_ls  = self._S_rot(T, self.ls)
        S_rot_hs  = self._S_rot(T, self.hs)

        S_spin_ls = self._S_spin_orbit(self.ls)
        S_spin_hs = self._S_spin_orbit(self.hs)

        S_conf_ls = self._S_config(self.ls)
        S_conf_hs = self._S_config(self.hs)

        dE_vib    = np.sum(E_vib_hs - E_vib_ls, axis=0)
        dS_vib    = np.sum(S_vib_hs - S_vib_ls, axis=0)

        dE_ele    = E_ele_hs  - E_ele_ls
        dS_ele    = S_ele_hs  - S_ele_ls
        dS_rot    = S_rot_hs  - S_rot_ls
        dS_spin   = S_spin_hs - S_spin_ls
        dS_conf   = S_conf_hs - S_conf_ls

        dU = dE_vib + dE_ele
        dS = dS_vib + dS_ele + dS_rot + dS_spin + dS_conf + k_times_n

        dH = dU + self.PdV
        dG = dH - dS*T

        return dG, dH, dS
#
### PRIVATE: Relative high-spin population for a given spin state
#
    def _nhs(self, spin_state:reader) -> float:

        nhs = np.count_nonzero(spin_state.spin != self.ls.spin)/self.centers

        if (nhs < 0.0) or (nhs > 1.0):
            raise RuntimeError(f"fractional high-spin population is {nhs:.2f}.")

        return nhs
#
### PRIVATE: Vibrational energy and entropy
#
    def _E_vib_and_S_vib(self, T:list, spin_state:reader, convert:float=1.0) -> tuple[np.ndarray[tuple[float], np.dtype[float]], np.ndarray[tuple[float], np.dtype[float]]]:

        frequencies = spin_state.frequencies

        T_size = T.size
        F_size = frequencies.size

        frq   = np.repeat(frequencies, T_size).reshape(F_size, T_size) 

        hkv   = h_over_k*frq

        hkvT  = hkv/T

        exp   = np.exp(-hkvT)

        oexp  = 1.0 - exp

        oexp  = np.where(oexp<log_thresh, log_thresh, oexp)

        loexp = np.log(oexp)

        E_vib = convert*k_times_n*(0.5*hkv + hkv*exp/oexp)
        S_vib = convert*k_times_n*(hkvT*exp/oexp - loexp)

        return E_vib, S_vib
#
### PRIVATE: Spin entropy
#
    def _S_spin_orbit(self, spin_state:reader, convert:float=1.0) -> float:

        idx     = self.hs.spin != self.ls.spin

        S_spin  = k_times_n*np.log(spin_state.jahnteller*(spin_state.spin[idx] + 1.0))

        S_orbit = k_times_n*np.log(2.0*spin_state.orbit[idx] + 1.0)

        S_soc   = np.sum(S_spin + S_orbit)*convert

        return S_soc
#
### PRIVATE: Entropy of mixing
#
    def _S_mix(self, nhs:np.ndarray[tuple[float], np.dtype[float]], d_nhs:int=0, convert:float=1.0) -> np.ndarray[tuple[float], np.dtype[float]]:

        if   d_nhs == 0:
            S_mix = np.where((nhs!=0.0)&(nhs!=1.0), -k_times_n*convert*(nhs*np.log(nhs) + (1.0 - nhs)*np.log(1.0 - nhs)), 0.0)

        elif d_nhs == 1:
            S_mix = np.where((nhs!=0.0)&(nhs!=1.0), k_times_n*convert*np.log((1.0 - nhs)/nhs), 0.0)

        elif d_nhs == 2:
            S_mix = np.where((nhs!=0.0)&(nhs!=1.0), k_times_n*convert/(nhs*(nhs - 1.0)), 0.0)

        else:
            raise RuntimeError(f"derivative order {d_nhs} not implemented for S_mix.")

        return S_mix
#
### PRIVATE: Entropy of configurations
#
    def _S_config(self, spin_state:reader) -> float:

        S_config = k_times_n*np.log(spin_state.configurations)

        return S_config
#
### PRIVATE: Electronic entropy
#
    def _S_ele(self, T:list, spin_state:reader, sigma:float=1e-02, convert:float=1.0) -> np.ndarray[tuple[float], np.dtype[float]]:

        orbital_energies = spin_state.orbital_energies

        S_ele = self._S_fermi_dirac(np.sort(orbital_energies), T, sigma=sigma, convert=convert)

        return S_ele
#
### PRIVATE: Rotational entropy
#
    def _S_rot(self, T:list, spin_state:reader) -> np.ndarray[tuple[float], np.dtype[float]]:

        rotational_symmetry    = spin_state.rotational_symmetry
        rotational_temperature = spin_state.rotational_temperature

        S_rot = np.zeros_like(T)

        if (not rotational_symmetry) or (not rotational_temperature): return S_rot

        T_cube = np.abs(T*T*T)
        theta  = rotational_symmetry*np.sqrt(rotational_temperature*T_cube)

        S_rot  = k_times_n*np.log(theta)

        return S_rot
#
### PRIVATE: Volume expansion
#
    def _PdV(self, convert:float=ang3_to_m3) -> float:

        PdV = atm_pressure*(self.hs.volume - self.ls.volume)*convert

        return PdV
#
### PRIVATE: Fermi-Dirac entropy
#
    @staticmethod
    def _S_fermi_dirac(energy:np.ndarray, T:list, sigma:float=1e-02, thresh:float=exp_thresh, convert:float=1.0) -> np.ndarray[tuple[float], np.dtype[float]]:

        T_size   = T.size
        ene_size = energy.size

        energy = np.repeat(energy, T_size).reshape(ene_size, T_size)
        fermi  = np.zeros_like(energy)
        dummy  = np.zeros_like(energy)

        dos = np.ones_like(energy) - sigma*energy**2
        dos = dos*(dos>0.0)

        E_over_KbT = energy/(k_times_n*T/ev_to_joule)

        fermi[E_over_KbT<-thresh] = 1.0
        fermi[E_over_KbT>-thresh] = 0.0

        idx = np.where( (E_over_KbT>-thresh) & (E_over_KbT<thresh) )
        fermi[idx] = 1.0/(1.0 + np.exp(E_over_KbT[idx]))

        idx = np.where( (fermi>0.0)&(fermi<1.0) )
        dummy[idx] = dos[idx]*(fermi[idx]*np.log(fermi[idx]) + (1.0 - fermi[idx])*np.log(1.0 - fermi[idx]))

        idx = np.logical_not(np.isnan(dummy))
        S_fermi_dirac = -k_times_n*np.trapz(energy[idx], dummy[idx], axis=0)*convert

        return S_fermi_dirac
#
### PRIVATE: Vibrational contribution to the heat capacity
#
    def _dE_vib_dT(self, T:list, spin_state:reader, convert:float=1.0) -> np.ndarray[tuple[float], np.dtype[float]]:

        frequencies = spin_state.frequencies

        T_size = T.size
        F_size = frequencies.size

        frq    = np.repeat(frequencies, T_size).reshape(F_size, T_size)

        hkv    = h_over_k*frq

        hkvT   = hkv/T

        exp    = np.exp(-hkvT)

        oexp   = 1.0 - exp

        hoexp  = hkvT/oexp

        dE_vib = convert*k_times_n*np.sum(exp*hoexp*hoexp, axis=0)

        return dE_vib
    
#############################################################################################################################

##     ## ####  ######  ########   #######   ######   ######   #######  ########  ####  ######  
###   ###  ##  ##    ## ##     ## ##     ## ##    ## ##    ## ##     ## ##     ##  ##  ##    ## 
#### ####  ##  ##       ##     ## ##     ## ##       ##       ##     ## ##     ##  ##  ##       
## ### ##  ##  ##       ########  ##     ##  ######  ##       ##     ## ########   ##  ##       
##     ##  ##  ##       ##   ##   ##     ##       ## ##       ##     ## ##         ##  ##       
##     ##  ##  ##    ## ##    ##  ##     ## ##    ## ##    ## ##     ## ##         ##  ##    ## 
##     ## ####  ######  ##     ##  #######   ######   ######   #######  ##        ####  ######

#############################################################################################################################

#############################################################################################################################