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
import multiprocessing as mp

from scipy.optimize import fsolve
#
### local libraries
#
from ..io.read         import reader
from .models           import Thermodynamic
from ..utils.constants import k_times_n
from ..utils.constants import spin_crossover_metals

cpus = mp.cpu_count()

#############################################################################################################################

########  ########  ######   ##     ## ##          ###    ########      ######   #######  ##       ##     ## ######## ####  #######  ##    ## 
##     ## ##       ##    ##  ##     ## ##         ## ##   ##     ##    ##    ## ##     ## ##       ##     ##    ##     ##  ##     ## ###   ## 
##     ## ##       ##        ##     ## ##        ##   ##  ##     ##    ##       ##     ## ##       ##     ##    ##     ##  ##     ## ####  ## 
########  ######   ##   #### ##     ## ##       ##     ## ########      ######  ##     ## ##       ##     ##    ##     ##  ##     ## ## ## ## 
##   ##   ##       ##    ##  ##     ## ##       ######### ##   ##            ## ##     ## ##       ##     ##    ##     ##  ##     ## ##  #### 
##    ##  ##       ##    ##  ##     ## ##       ##     ## ##    ##     ##    ## ##     ## ##       ##     ##    ##     ##  ##     ## ##   ### 
##     ## ########  ######    #######  ######## ##     ## ##     ##     ######   #######  ########  #######     ##    ####  #######  ##    ##

#############################################################################################################################

class RegularSolution(Thermodynamic):
    """
    Thermodynamic regular solution model.

    .. [1] J. H. Hildebrand and R. L. Scott, Regular Solutions; Prentice-Hall International Series in Chemistry; Prentice-Hall, 1962.
    .. [2] J. H. Hildebrand, J. M. Prausnitz, and R. L. Scott, Regular and Related Solutions: the Solubility of Gases, Liquids,
           and Solids; Van Nostrand Reinhold Co.: New York, 1970.

    Parameters
    ----------
        ls : reader
            Defines the reference low-spin state read from an electronic structure output.

        hs : reader
            Defines the reference high-spin state read from an electronic structure output.

        centers : float, default=1.0
            Number of spin-crossover metal centers used to normalize reported data.

        metals : tuple(str), default=('Cr', 'Mn', 'Fe', 'Co')
            Define the transition metals of interest for automatic determination of the parameter 'centers'.

    Modules
    -------
        spin_crossover_energy
            Compute the spin-crossover energy, or spin gap.

        transition_temperature
            Compute the transition temperature that corresponds to a high-spin population of 1/2.

        high_spin_population
            Compute the evolution of the the relative high-spin population, or magnetization.

        gibbs_free_energy
            Compute the Gibbs free energy for a given isotherm.

    Model Modules
    -------------
        heat_capacity
            Compute the thermal evolution of the heat capacity.

    Fine Control
    ------------
        set_nhs_list
            Override the default list of points with a custom module.
        
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
#
### Initialize
#
    def __init__(self, ls:reader, hs:reader, centers:float=0.0, metals=spin_crossover_metals) -> None:
        super().__init__(ls, hs, centers=centers, metals=metals)
#
### PUBLIC: Compute the spin-crossover energy
#
    def spin_crossover_energy(self, zero_point_energy:bool=False) -> float:
        """
        Compute the spin-crossover energy between a high- and the low-spin state.

        Parameters
        ----------
            zero_point_energy : bool, default=False
                Define whether to include the zero-point energy correction.

        Returns
        -------
            delta_energy : float
                Spin-crossover energy normalized to the number of spin transition centers.

        Units
        -----
        kJ/mol
        """

        delta_energy = 1e-03*self.spin_gap(zero_point_energy=zero_point_energy) / self.centers

        return delta_energy
#
### PUBLIC: Compute the transition temperature
#
    def transition_temperature(self, guess:float=298.15, tol:float=1e-03) -> tuple[float, float, float]:
        """
        Compute the transition temperature for the spin conversion.

        Parameters
        ----------
            guess : float, default=298.15 K
                Define the initial guess for the numeric procedure.

            tol : float, default=1.0E-03 K
                Define the accuracy threshold for the transition temperature.

        Returns
        -------
            T_half : float
                Transition temperature.

            dH : float
                Enthalpy change.

            dS: float
                Entropy change.

        Units
        -----
            K, kJ/mol, J/(mol K), respectively
        """

        T_half, dH, dS = self._temperature(0.5, guess=guess, tol=tol)

        return T_half.item(), dH.item(), dS.item()
#
### PUBLIC: Compute relative high-spin population as function of T
#
    def high_spin_population(self, points:int=128, guess:float=298.15, tol:float=1e-03) -> dict[str, np.ndarray[tuple[float], np.dtype[float]]]:
        """
        Compute the thermal evolution of the relative high-spin population, or magnetization.

        ΔH = T ΔS + k_B N_A T ln[(1 - n_HS)/n_HS]

        Parameters
        ----------
            points : int, default=128
                Specify the total number of points for the curve.

            guess : float, default=298.15 K
                Define the initial guess for the numeric procedure.

            tol : float, default=1.0E-03 K
                Define the accuracy threshold for the associated temperature.

        Returns
        -------
            'nhs' : np.ndarray[ tuple[float], np.dtype[float] ]
                Relative high-spin population.

            'T' : np.ndarray[ tuple[float], np.dtype[float] ]
                Temperature.

            'dH' : np.ndarray[ tuple[float], np.dtype[float] ]
                Enthalpy change.

            'dS' : np.ndarray[ tuple[float], np.dtype[float] ]
                Entropy change.

            'dG' : np.ndarray[ tuple[float], np.dtype[float] ]
                Gibbs free energy change.

        Units
        -----
            dimensionless, K, kJ/mol, J/(mol K), kJ/mol, respectively
        """

        nhs = self.set_nhs_list(points, self.min_nhs, self.max_nhs)

        pool  = mp.Pool(cpus)
        mparg = [ (nhs, guess, tol) for nhs in nhs ]

        T_dH_dS = pool.starmap_async(self._temperature, mparg).get()
        pool.close()

        T_dH_dS = np.asarray(T_dH_dS).squeeze(2)

        keepidx = (T_dH_dS[:,0] != guess)

        T_dH_dS = T_dH_dS[keepidx]

        dG, _,_ = self.delta_gibbs(T_dH_dS[:,0])

        nhs_T_dH_dS_dG = {
            "nhs" : nhs[keepidx],
            "T"   : T_dH_dS[:,0],
            "dH"  : T_dH_dS[:,1],
            "dS"  : T_dH_dS[:,2],
            "dG"  : 1e-03*dG
        }

        return nhs_T_dH_dS_dG
#
### PUBLIC: Compute Gibbs free energy as function of T (isothermal or computed)
#
    def gibbs_free_energy(self, temperature:float, points:int=128) -> dict[str, np.ndarray[tuple[float], np.dtype[float]]]:
        """
        Compute the Gibbs free energy for a given isotherm.

        G = n_HS G_HS + (1 - n_HS) G_LS + T S_mix

        Parameters
        ----------           
            temperature : float
                Define the temperature for the isotherm.

            points : int, default=128
                Specify the total number of points for the magnetization curve.

        Returns
        -------
            'nhs' : np.ndarray[ tuple[float], np.dtype[float] ]
                Relative high-spin population.

            'G' : np.ndarray[ tuple[float], np.dtype[float] ]
                Gibbs free energy.

        Units
        -----
        dimensionless, kJ/mol, respectively
        """
        
        if temperature <= 0.0:
            raise RuntimeError("temperature must be positive and larger than zero.")

        nhs = self.set_nhs_list(points, self.min_nhs, self.max_nhs)
        T   = np.repeat(temperature, nhs.size)

        dG,_,_ = self.delta_gibbs(T)

        Gibbs  = nhs*dG - self._S_mix(nhs)*T

        nhs_G  = {
            "nhs" : nhs,
            "G"   : 1e-03*Gibbs
        }

        return nhs_G
#
### PUBLIC: Compute the heat capacity as function of temperature
#
    def heat_capacity(self, points:int=128, guess:float=298.15, tol:float=1e-03) -> dict[str, np.ndarray[tuple[float], np.dtype[float]]]:
        """
        Compute the thermal evolution of the heat capacity.

        Parameters
        ----------           
            points : int, default=128
                Specify the total number of points for the curve.
            
            guess : float, default=298.15 K
                Define the initial guess for the numeric procedure.

            tol : float, default=1.0E-03 K
                Define the accuracy threshold for the associated temperature.

        Returns
        -------
            'T' : np.ndarray[ tuple[float], np.dtype[float] ]
                Temperature.

            'Cp' : np.ndarray[ tuple[float], np.dtype[float] ]
                Heat capacity.

        Units
        -----
        K, J/(mol K), respectively
        """

        nhs_points = self.set_nhs_list(points, self.min_nhs, self.max_nhs)

        pool  = mp.Pool(cpus)
        mparg = [ (nhs, guess, tol) for nhs in nhs_points ]

        T_dH_dS = pool.starmap_async(self._temperature, mparg).get()
        pool.close()

        T_dH_dS = np.asarray(T_dH_dS).squeeze(2)

        keepidx = (T_dH_dS[:,0] != guess)

        nhs     = nhs_points[keepidx]
        T_dH_dS = T_dH_dS[keepidx]

        T, dH   = T_dH_dS[:,0], T_dH_dS[:,1]
        dG,_,_  = self.delta_gibbs(T)

        dE_vib_ls = self._dE_vib_dT(T, self.ls)
        dE_vib_hs = self._dE_vib_dT(T, self.hs)

        dGkT  = dG/(k_times_n*T)
        kTT   = 1.0/(k_times_n*T*T)

        exp   = np.exp(-dGkT)
        oexp  = 1.0/(1.0 + exp)

        Cp    = (1.0 - nhs)*dE_vib_ls + nhs*dE_vib_hs + kTT*dH*dH*exp*oexp*oexp

        T_Cv  = {
            "T"  : T,
            "Cp" : Cp,
        }

        return T_Cv
#
### PRIVATE: Compute transition temperature
#
    def _temperature(self, nhs, guess:float=298.15, tol:float=1e-03) -> tuple[np.ndarray[tuple[float], np.dtype[float]], np.ndarray[tuple[float], np.dtype[float]], np.ndarray[tuple[float], np.dtype[float]]]:
        """
        Compute the temperature associated to a specific value of the relative high-spin population, or magnetization, n_HS.

        ΔH = T\*ΔS + k_B\*N_A\*T\*ln[(1 - n_HS)/n_HS]
        """

        if isinstance(nhs, float) and isinstance(guess, float):
            nhs   = np.asarray([nhs])
            guess = np.asarray([guess])

        T_half = fsolve(self._T_from_nhs, guess, args=nhs, xtol=tol)

        _, dH, dS = self.delta_gibbs(T_half)

        return T_half, 1e-03*dH, dS
#
### PRIVATE: Determine temperature fron high-spin population
#
    def _T_from_nhs(self, T:list, *args:tuple) -> np.ndarray[tuple[float], np.dtype[float]]:
        """
        Find the extremum T given n_HS for the equation

        ΔG - k_B N_A T ln[(1 - n_HS)/n_HS] = 0
        """

        dG,_,_ = self.delta_gibbs(T)
        S_mix  = self._S_mix(*args, d_nhs=1)

        Gibbs  = dG - S_mix*T

        return Gibbs
    
#############################################################################################################################

 ######  ##       ####  ######  ##     ## ######## ######## ########     ########  ########  ####  ######  ##    ##    ###    ##     ## ######## ########  
##    ## ##        ##  ##    ## ##     ##    ##    ##       ##     ##    ##     ## ##     ##  ##  ##    ## ##   ##    ## ##   ###   ### ##       ##     ## 
##       ##        ##  ##       ##     ##    ##    ##       ##     ##    ##     ## ##     ##  ##  ##       ##  ##    ##   ##  #### #### ##       ##     ## 
 ######  ##        ##  ##       #########    ##    ######   ########     ##     ## ########   ##  ##       #####    ##     ## ## ### ## ######   ########  
      ## ##        ##  ##       ##     ##    ##    ##       ##   ##      ##     ## ##   ##    ##  ##       ##  ##   ######### ##     ## ##       ##   ##   
##    ## ##        ##  ##    ## ##     ##    ##    ##       ##    ##     ##     ## ##    ##   ##  ##    ## ##   ##  ##     ## ##     ## ##       ##    ##  
 ######  ######## ####  ######  ##     ##    ##    ######## ##     ##    ########  ##     ## ####  ######  ##    ## ##     ## ##     ## ######## ##     ##

#############################################################################################################################

class SlichterDrickamer(Thermodynamic):
    """
    Thermodynamic extension to the regular soultion model proposed by Slichter and Drickamer.

    .. [1] C. P. Slichter and H. G. Drickamer J. Chem. Phys. 56, 2142-2160 (1972).

    Parameters
    ----------
        ls : reader
            Defines the reference low-spin state read from an electronic structure output.

        hs : reader
            Defines the reference high-spin state read from an electronic structure output.

        centers : float, default=1.0
            Number of spin-crossover metal centers used to normalize reported data.

        metals : tuple(str), default=('Cr', 'Mn', 'Fe', 'Co')
            Define the transition metals of interest for automatic determination of the parameter 'centers'.

    Modules
    -------
        spin_crossover_energy
            Compute the spin-crossover energy, or spin gap.

        transition_temperature
            Compute the transition temperature that corresponds to a high-spin population of 1/2.

        high_spin_population
            Compute the evolution of the the relative high-spin population, or magnetization.

        gibbs_free_energy
            Compute the Gibbs free energy for a given isotherm.

    Model Modules
    -------------
        fit_boltzmann_interaction
            Fit the interaction parameter using Boltzmann weights for a series of intermediate spin configurations.

        fit_leastsqrs_interaction
            Fit the interaction parameter from a least squares fit for a series intermediate spin configurations.

    Fine Control
    ------------
        set_nhs_list
            Override the default list of points with a custom module.
        
        set_min_nhs
            Set the lower threshold for the relative high-spin population.

        set_max_nhs
            Set the upper threshold for the relative high-spin population.

        set_interaction
            Set the default value of the interaction parameter for use during runtime.

    Attributes
    ----------
        min_nhs : float, default=1.00E-06
            Lower bound for the relative high-spin population of interest.

        max_nhs : float, default=9.99E-01
            Upper bound for the relative high-spin population of interest.

        interaction : float, default=0.0
            Phenomenological interaction parameter.
    """
#
### Initialize
#
    def __init__(self, ls:reader, hs:reader, interaction:float=0.0, centers:float=0.0, metals=spin_crossover_metals) -> None:
        super().__init__(ls, hs, centers=centers, metals=metals)

        self.interaction = 1e+03*interaction
#
### PUBLIC: Compute the spin-crossover energy
#
    def spin_crossover_energy(self, zero_point_energy:bool=False) -> float:
        """
        Compute the spin-crossover energy between a high- and the low-spin state.

        Parameters
        ----------
            zero_point_energy : bool, default=False
                Define whether to include the zero-point energy correction.

        Returns
        -------
            delta_energy : float
                Spin-crossover energy normalized to the number of spin transition centers.

        Units
        -----
        kJ/mol
        """

        delta_energy = 1e-03*self.spin_gap(zero_point_energy=zero_point_energy) / self.centers

        return delta_energy
#
### PUBLIC: Compute the transition temperature
#
    def transition_temperature(self, guess:float=298.15, tol:float=1e-03, interaction:float=0.0) -> tuple[float, float, float]:
        """
        Compute the transition temperature for the spin conversion.

        Parameters
        ----------
            guess : float, default=298.15 K
                Define the initial guess for the numeric procedure.

            tol : float, default=1.0E-03 K
                Define the accuracy threshold for the transition temperature.

            interaction : float, default=0.0 kJ/mol
                Define the phenomenological interaction parameter.
                When provided, it is given preference upon the default value, should it exist.

        Returns
        -------
            T_half : float
                Transition temperature.

            dH : float
                Enthalpy change.

            dS: float
                Entropy change.

        Units
        -----
            K, kJ/mol, J/(mol K), respectively
        """

        interaction    = self.interaction if not interaction else 1e+03*interaction

        T_half, dH, dS = self._temperature(0.5, guess=guess, tol=tol, interaction=interaction)

        return T_half.item(), dH.item(), dS.item()
#
### PUBLIC: Compute relative high-spin population as function of T
#
    def high_spin_population(self, points:int=128, guess:float=298.15, tol:float=1e-03, interaction:float=0.0) -> dict[str, np.ndarray[tuple[float], np.dtype[float]]]:
        """
        Compute the thermal evolution of the relative high-spin population, or magnetization, n_HS.

        ΔH + Γ (1 - 2 n_HS) = T ΔS + k_B N_A T ln[(1 - n_HS)/n_HS]

        Parameters
        ----------
            points : int, default=128
                Specify the total number of points for the curve.

            guess : float, default=298.15 K
                Define the initial guess for the numeric procedure.

            tol : float, default=1.0E-03 K
                Define the accuracy threshold for the associated temperature.

            interaction : float, default=0.0 kJ/mol
                Define the phenomenological interaction parameter.
                When provided, it is given preference upon the default value, should it exist.

        Returns
        -------
            'nhs' : np.ndarray[ tuple[float], np.dtype[float] ]
                Relative high-spin population.

            'T' : np.ndarray[ tuple[float], np.dtype[float] ]
                Temperature.

            'dH' : np.ndarray[ tuple[float], np.dtype[float] ]
                Enthalpy change.

            'dS' : np.ndarray[ tuple[float], np.dtype[float] ]
                Entropy change.

            'dG' : np.ndarray[ tuple[float], np.dtype[float] ]
                Gibbs free energy change.

        Units
        -----
            dimensionless, K, kJ/mol, J/(mol K), kJ/mol, respectively
        """

        interaction = self.interaction if not interaction else 1e+03*interaction

        nhs_points  = self.set_nhs_list(points, self.min_nhs, self.max_nhs)

        pool  = mp.Pool(cpus)
        mparg = [ (nhs, guess, tol, interaction) for nhs in nhs_points ]

        T_dH_dS = pool.starmap_async(self._temperature, mparg).get()
        pool.close()

        T_dH_dS = np.asarray(T_dH_dS).squeeze(2)

        keepidx = (T_dH_dS[:,0] != guess)

        T_dH_dS = T_dH_dS[keepidx]

        dG, _,_ = self.delta_gibbs(T_dH_dS[:,0])

        nhs_T_dH_dS_dG = {
            "nhs" : nhs_points[keepidx],
            "T"   : T_dH_dS[:,0],
            "dH"  : T_dH_dS[:,1],
            "dS"  : T_dH_dS[:,2],
            "dG"  : 1e-03*dG
        }

        return nhs_T_dH_dS_dG
#
### PUBLIC: Compute Gibbs free energy as function of T (isothermal or computed)
#
    def gibbs_free_energy(self, temperature:float, points:int=128, interaction:float=0.0) -> dict[str, np.ndarray[tuple[float], np.dtype[float]]]:
        """
        Compute the Gibbs free energy for a given isotherm.

        G = n_HS G_HS + (1 - n_HS) G_LS + Γ n_HS (1 - n_HS) + T S_mix

        Parameters
        ----------           
            temperature : float
                Define the temperature for the isotherm

            points : int, default=128
                Specify the total number of points for the curve.

            interaction : float, default=0.0 kJ/mol
                Define the phenomenological interaction parameter.
                When provided, it is given preference upon the default value, should it exist.

        Returns
        -------
            'nhs' : np.ndarray[ tuple[float], np.dtype[float] ]
                Relative high-spin population.

            'G' : np.ndarray[ tuple[float], np.dtype[float] ]
                Gibbs free energy.

        Units
        -----
        dimensionless, kJ/mol, respectively
        """

        interaction = self.interaction if not interaction else 1e+03*interaction
        
        if temperature <= 0.0:
            raise RuntimeError("temperature must be positive and larger than zero.")
        
        nhs = self.set_nhs_list(points, self.min_nhs, self.max_nhs)
        T   = np.repeat(temperature, nhs.size)

        dG,_,_ = self.delta_gibbs(T)

        Gibbs  = nhs*dG + interaction*nhs*(1.0 - nhs) - self._S_mix(nhs)*T

        nhs_G  = {
            "nhs" : nhs,
            "G"   : 1e-03*Gibbs
        }

        return nhs_G
#
### PUBLIC: Set default interaction parameter
#
    def set_interaction(self, interaction:float) -> None:
        """
        Fix the default value of the interaction parameter for use during runtime.

        Parameters
        ----------
            interaction : float, kJ/mol
                Phenomenological interaction parameter.

        Returns
        -------
            None
        """

        interaction *= 1e-03

        self.interaction = interaction
#
### PUBLIC: Fit interaction parameter using Botlzmann weights
#
    def fit_boltzmann_interaction(self, ms:list[reader], guess:float=298.15, tol:float=1e-03, interaction:float=0.0) -> tuple[float, list[float]]:
        """
        Fit the interaction parameter from Boltzmann weights using a set of intermediate spin configurations.

        Parameters
        ----------
            ms : list[reader]
                Define a list of intermediate spin configurations.

            guess : float, default=298.15 K
                Define the initial guess for the numeric procedure.

            tol : float, default=1.0E-03 K
                Define the accuracy threshold for the associated temperature.

            interaction : float, default=0.0 kJ/mol
                Define the phenomenological interaction parameter.
                When provided, it is given preference upon the default value, should it exist.

        Returns
        -------
            interaction : float
                Phenomenological interaction parameter.

            weights : list[float]
                Associated weight for each intermediate spin configuration.

        Units
        -----
            kJ/mol, dimensionless, respectively
        """

        interaction = self.interaction if not interaction else 1e+03*interaction
        
        states = len(ms)
        ncpus  = states if cpus > states else cpus

        pool   = mp.Pool(ncpus)
        mparg  = [(spin_state, guess, tol, interaction) for spin_state in ms]

        Eml_nhs_factor_c = pool.starmap_async(self._intermediate_spin_data, mparg).get()
        pool.close()

        Eml_nhs_factor_c = np.asarray(Eml_nhs_factor_c)

        fc = np.asarray(Eml_nhs_factor_c[:,2])
        gc = np.asarray(Eml_nhs_factor_c[:,3])

        weights = fc/np.sum(fc)

        interaction = np.sum(weights*gc)

        return 1e-03*interaction, weights.tolist()
#
### PUBLIC: Fit interaction parameter using least squares
#
    def fit_leastsqrs_interaction(self, ms:list[reader], guess:float=298.15, tol:float=1e-03, interaction:float=0.0) -> tuple[float, float]:
        """
        Fit the interaction parameter from least squares using a set of intermediate spin configurations.

        Parameters
        ----------
            ms : list[reader]
                Define a list of intermediate spin configurations.

            guess : float, default=298.15 K
                Define the initial guess for the numeric procedure.

            tol : float, default=1.0E-03 K
                Define the accuracy threshold for the associated temperature.

            interaction : float, default=0.0 kJ/mol
                Define the phenomenological interaction parameter.
                When provided, it is given preference upon the default value, should it exist.

        Returns
        -------
            interaction : float
                Phenomenological interaction parameter.

            R2 : float
                Coefficient of determination.

        Units
        -----
            kJ/mol, dimensionless, respectively
        """

        interaction = self.interaction if not interaction else 1e+03*interaction
        
        states = len(ms)
        ncpus  = states if cpus > states else cpus

        pool   = mp.Pool(ncpus)
        mparg  = [(spin_state, guess, tol, interaction) for spin_state in ms]

        Eml_nhs_factor_c = pool.starmap_async(self._intermediate_spin_data, mparg).get()
        pool.close()

        Eml_nhs_factor_c = np.asarray(Eml_nhs_factor_c)

        E_hl = self.spin_gap(zero_point_energy=True)
        E_ml = np.asarray(Eml_nhs_factor_c[:,0])
        nhs  = np.asarray(Eml_nhs_factor_c[:,1])

        onenhs  = nhs*(1.0 - nhs)
        nhsE_hl = nhs*E_hl

        interaction = np.sum((E_ml - nhsE_hl)*onenhs)/np.sum(onenhs*onenhs)

        fit = nhsE_hl + onenhs*interaction
        st2 = E_ml - np.sum(E_ml)/E_ml.size
        se2 = E_ml - fit

        sst = st2*st2
        sse = se2*se2

        R2  = 1.0 - np.sum(sse)/np.sum(sst)

        return 1e-03*interaction, R2
#
### PRIVATE: Compute transition temperature
#
    def _temperature(self, nhs, guess:float=298.15, tol:float=1e-03, interaction:float=0.0) -> tuple[np.ndarray[tuple[float], np.dtype[float]], np.ndarray[tuple[float], np.dtype[float]], np.ndarray[tuple[float], np.dtype[float]]]:
        """
        Compute the temperature associated to a specific value of the relative high-spin population, or magnetization, n_HS.

        ΔH + Γ\*(1 - 2\*n_HS) = T\*ΔS + k_B\*N_A\*T\*ln[(1 - n_HS)/n_HS]
        """

        interaction = self.interaction if not interaction else 1e-03*interaction

        if isinstance(nhs, float) and isinstance(guess, float):
            nhs   = np.asarray([nhs])
            guess = np.asarray([guess])

        T_half = fsolve(self._T_from_nhs, guess, args=(nhs, interaction), xtol=tol)

        _, dH, dS = self.delta_gibbs(T_half)

        return T_half, 1e-03*dH, dS
#
### PRIVATE: Determine temperature fron high-spin population
#
    def _T_from_nhs(self, T:list, *args:tuple) -> np.ndarray[tuple[float], np.dtype[float]]:
        """
        Find the extremum T given n_HS for the equation

        ΔG + Γ\*(1 - 2\*n_HS) - k_B\*N_A\*T\*ln[(1 - n_HS)/n_HS] = 0
        """

        nhs, interaction = args

        interaction = self.interaction if not interaction else 1e+03*interaction

        dG,_,_ = self.delta_gibbs(T)
        S_mix  = self._S_mix(nhs, d_nhs=1)

        Gibbs  = dG + interaction*(1.0 - 2.0*nhs) - S_mix*T

        return Gibbs
#
### PRIVATE: Calculate diverse data for an intermediate spin state
#
    def _intermediate_spin_data(self, spin_state:reader, guess:float=298.15, tol:float=1e-03, interaction:float=0.0) -> tuple[float, np.ndarray[tuple[float], np.dtype[float]], np.ndarray[tuple[float], np.dtype[float]], np.ndarray[tuple[float], np.dtype[float]]]:

        interaction = self.interaction if not interaction else 1e+03*interaction

        nhs   = self._nhs(spin_state)

        E_hl  = self.spin_gap(zero_point_energy=True)
        E_ml  = self.spin_gap(spin_state, zero_point_energy=True)

        dE_ml = E_ml - nhs*E_hl

        T,_,_ = self._temperature(nhs, guess=guess, tol=tol, interaction=interaction)

        KbT   = k_times_n*T

        boltzmann_factor = np.exp(-dE_ml/KbT)
        interaction_c    = dE_ml/(nhs*(1.0 - nhs))
        
        return E_ml, nhs, boltzmann_factor.item(), interaction_c

#############################################################################################################################
