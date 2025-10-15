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
#
### local libraries
#
from .constants import spin_crossover_metals
from .constants import ev_to_joule, ang3_to_m3
from .constants import k_times_n, atm_pressure

from .utils     import extract_data
from .utils     import calculate_thalf, calculate_nhs, calculate_interaction, calculate_rsm

from .states    import generate_mixed_spin_states

#############################################################################################################################

 ######  ##          ###     ######   ######  ########  ######
##    ## ##         ## ##   ##    ## ##    ## ##       ##    ##
##       ##        ##   ##  ##       ##       ##       ##
##       ##       ##     ##  ######   ######  ######    ######
##       ##       #########       ##       ## ##             ##
##    ## ##       ##     ## ##    ## ##    ## ##       ##    ##
 ######  ######## ##     ##  ######   ######  ########  ######

#############################################################################################################################
### Define specific spin state

class SpinState:

    __slots__ = ("states", "crossovers", "nhs", "step")
#
## Initialize
#   
    def __init__(self, states:int=2) -> None:
        
        self.states     = states
        self.crossovers = self.states - 1
        
        self.nhs    = []
        self.step   = []
#
### Prepare tuple for use as args
#
    def zipit(self, idx=0, interaction:float=0.0, nhs:float=None) -> tuple:

        steps = len(self.step)

        if not -1 < idx < steps: raise RuntimeError(f"index {idx} out of range [0-{steps-1}]")
            
        if not self.step[idx]: raise RuntimeError(f"spin-step with index {idx} is empty")
            
        crossover_energy, spin_entropy, volume_expansion, frequencies, orbitalenergies, rotational = self.step[idx]

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
        
#############################################################################################################################
### Thermodynamic properties

class ThermodynamicProperties:
#
### Initialize
#
    def __init__(self, centers:float=0.0) -> None:
        
        self.centers = centers
        
        self.atomlist         = None
        self.ls_magnetization = None
        self.hs_magnetization = None
        self.ls_atoms         = None
        self.hs_atoms         = None
        self.ls_elements      = None
        self.hs_elements      = None
        self.ls_jahnteller    = None
        self.hs_jahhteller    = None
        self.ls_angular       = None
        self.hs_angular       = None
#
### Set spin-state reference data
#
    def setref(self, ls:any=None, hs:any=None, metals:tuple=spin_crossover_metals) -> None:
        
        if (ls is None) and (hs is None):
            raise RuntimeError(f"low- nor high-spin reference selected")
        
        if ls is not None:
            self.ls_magnetization = ls.magnetization
            self.ls_atoms         = ls.atoms
            self.ls_elements      = ls.elements
            
            self.ls_jahnteller = ls.jahnteller
            self.ls_angular    = np.array(ls.magnetization.size*[ls.orbit])
            
        if hs is not None:
            self.hs_magnetization = hs.magnetization
            self.hs_atoms         = hs.atoms
            self.hs_elements      = hs.elements
            
            self.hs_jahnteller = hs.jahnteller
            self.hs_angular    = np.array(hs.magnetization.size*[hs.orbit])
            
        if (ls is not None) and (hs is not None):
            self.atomlist = np.repeat(ls.elements, ls.atoms)
            
            if not self.centers:
                self.centers = np.count_nonzero(np.isin(self.atomlist[ls.magnetization != hs.magnetization], metals))
#
### Relative high-spin population
#
    def nhs(self, ls:any=None, hs:any=None) -> float:
        
        ls_frac = np.count_nonzero(ls.magnetization != self.ls_magnetization)
        hs_frac = np.count_nonzero(hs.magnetization != self.ls_magnetization)
        
        fraction = (hs_frac - ls_frac)/self.centers
        
        if (fraction < 0.0) or (fraction > 1.0):
            raise RuntimeError(f"fractional high-spin population is {fraction:.2f}")
        
        return fraction
#
### Spin-crossover energy
#
    def sco(self, ls=None, hs=None, zpe=True, convert:float=ev_to_joule) -> float:

        energy = (hs.energy - ls.energy)*convert

        if zpe: energy += (hs.zeropointenergy - ls.zeropointenergy)*convert
        
        return energy
#
### Spin entropy
#
    def s_spin(self, ls:any=None, hs:any=None, convert:float=1.0) -> float:
        
        lsdiff = ls.magnetization == self.ls_magnetization
        hsdiff = hs.magnetization == self.ls_magnetization
        idx    = lsdiff != hsdiff
        
        spin   = k_times_n*(np.log((hs.magnetization[idx] + 1.0)*self.hs_jahnteller) - np.log((ls.magnetization[idx] + 1.0)*self.ls_jahnteller))
        
        orbit  = k_times_n*(np.log(2.0*self.hs_angular[idx] + 1.0) - np.log(2.0*self.ls_angular[idx] + 1.0))

        spin_orbit = np.sum(spin + orbit)*convert

        return spin_orbit
#
### Entropy of mixing
#
    def s_mix(self, nhs:float=0.0, convert:float=1.0) -> float:

        mixing = -k_times_n*convert*(nhs*np.log(nhs) + (1.0 - nhs)*np.log(1.0 - nhs)) if not nhs == 1.0 else 0.0
                        
        return mixing
#
### Volume expansion
#
    def pdv(self, ls:any=None, hs:any=None, convert:float=atm_pressure*ang3_to_m3) -> float:

        expansion = (hs.volume - ls.volume)*convert

        return expansion

#############################################################################################################################
### Phenomenological interaction parameter

class PhenomenologicalInteractionParameter:
#
### Initialize
#
    def __init__(self, e_hl:float=0.0, e_ml:np.ndarray=None, nhs:np.ndarray=None, boltzmann_factor:np.ndarray=None, interaction_set:np.ndarray=None) -> None:
        
        self.E_hl             = e_hl
        self.E_ml             = np.asarray(e_ml, dtype=float)
        self.nhs              = np.asarray(nhs, dtype=float)
        self.interaction_set  = np.asarray(interaction_set, dtype=float)
        self.boltzmann_factor = np.asarray(boltzmann_factor, dtype=float)
        
        self.interaction_boltzmann  = 0.0
        self.interaction_leastsqrs  = 0.0
#
### Interaction parameter using Boltzmann weights
#
    def boltzmann(self, weights:bool=False) -> float:
        
        fc    = self.boltzmann_factor
        gc    = self.interaction_set
        
        wc    = fc/np.sum(fc)
        
        interaction = np.sum(wc*gc)
        
        self.interaction_boltzmann = interaction
        
        if weights: return interaction, wc

        return interaction
#
### Interaction parameter using least squares
#
    def leastsqrs(self, detcoef:bool=False) -> float:
        
        E_hl  = self.E_hl
        E_ml  = self.E_ml
        nhs   = self.nhs
        
        onehs = nhs*(1.0 - nhs)
        ehlhs = nhs*E_hl
        
        interaction = np.sum((E_ml - ehlhs)*onehs)/np.sum(onehs*onehs)
        
        if detcoef:
            fit = ehlhs + interaction*onehs
            st2 = E_ml - np.sum(E_ml)/E_ml.size
            se2 = E_ml - fit
            
            sst = st2*st2
            sse = se2*se2
            
            r2  = 1.0 - np.sum(sse)/np.sum(sst)
            
        self.interaction_leastsqrs = interaction

        if detcoef: return interaction, r2
        
        return interaction

#############################################################################################################################

##     ##  #######  ########  ##     ## ##       ########  ######
###   ### ##     ## ##     ## ##     ## ##       ##       ##    ##
#### #### ##     ## ##     ## ##     ## ##       ##       ##
## ### ## ##     ## ##     ## ##     ## ##       ######    ######
##     ## ##     ## ##     ## ##     ## ##       ##             ##
##     ## ##     ## ##     ## ##     ## ##       ##       ##    ##
##     ##  #######  ########   #######  ######## ########  ######

#############################################################################################################################
### Calculate spin-crossover energy

def spin_crossover_energy(ls:any=None, hs:any=None, centers:float=0.0, metals:tuple=spin_crossover_metals, zero_point_energy=False) -> float:

    Thermo = ThermodynamicProperties(centers=centers)

    Thermo.setref(ls=ls, hs=hs, metals=metals)

    transition_energy = Thermo.sco(ls=ls, hs=hs, zpe=zero_point_energy)

    return 1e-3*transition_energy/Thermo.centers

#############################################################################################################################
### Calculate transition temperature

def transition_temperature(ls:any=None, hs:any=None, centers:float=0.0, metals:tuple=spin_crossover_metals, dH_and_dS=False, guess=298.15) -> float:

    States = SpinState(states=2)
    Thermo = ThermodynamicProperties(centers=centers)

    Thermo.setref(ls=ls, hs=hs, metals=metals)

    step = extract_data(ls=ls, hs=hs, thermo=Thermo)

    States.step += [step]
    
    transition_temperature, delta_enthalpy, delta_entropy = calculate_thalf(data=States, guess=guess, centers=Thermo.centers)

    if not dH_and_dS:
        return transition_temperature
    
    return transition_temperature, 1e-03*delta_enthalpy, delta_entropy
    
#############################################################################################################################
### Calculate relative high-spin population

def high_spin_population(ls:any=None, hs:any=None, points:list|int=256, interaction:float=0.0, centers:float=0.0, metals:tuple=spin_crossover_metals) -> np.ndarray:

    if interaction: interaction *= 1e3

    States = SpinState(states=2)
    Thermo = ThermodynamicProperties(centers=centers)

    Thermo.setref(ls=ls, hs=hs, metals=metals)

    step, nhs = extract_data(ls=ls, hs=hs, thermo=Thermo, return_nhs=True)

    States.step += [step]
    States.nhs  += [nhs]

    option = "direct" if interaction <= 0.0 else "numeric"

    temperature_and_nhs = calculate_nhs(option=option, guess=298.15, points=points, data=States, interaction=interaction)

    return temperature_and_nhs

#############################################################################################################################
### Calculate the phenomenological interaction parameter

def interaction_parameter(ls:any=None, hs:any=None, ms:list=None, interaction:float=0.0, centers:float=0.0,
                          boltzmann:bool=False, least_squares:bool=True, option:str="scf",
                          metals:tuple=spin_crossover_metals) -> tuple[float, float]:
    
    if interaction: interaction *= 1e3
    
    if isinstance(ms,str):
        ms = generate_mixed_spin_states(ls=ls, hs=hs, ms=ms)

    States = SpinState(states=2+len(ms))
    Thermo = ThermodynamicProperties(centers=centers)

    Thermo.setref(ls=ls, hs=hs, metals=metals)

    step = extract_data(ls=ls, hs=hs, thermo=Thermo)

    States.step += [step]

    configurations = calculate_interaction(ls=ls, hs=hs, ms=ms, interaction=interaction, states=States, thermo=Thermo, option=option)

    SlichterDrickamerModel = PhenomenologicalInteractionParameter(e_hl=step[0],
                                                                  e_ml=configurations[:,0],
                                                                  nhs=configurations[:,1],
                                                                  boltzmann_factor=configurations[:,2],
                                                                  interaction_set=configurations[:,3])
    if least_squares: boltzmann = False

    if boltzmann:
        interaction, weights = SlichterDrickamerModel.boltzmann(weights=True)

        return 1e-3*interaction, weights
    
    interaction, determination_coefficient = SlichterDrickamerModel.leastsqrs(detcoef=True)

    return 1e-3*interaction, determination_coefficient

#############################################################################################################################
### Calculate Gibbs as a function of high-spin population

def regular_solution_model(ls:any=None, hs:any=None, ms:list=None, temperature:float=0.0, points:list|int=256,
                           interaction:float=0.0, centers:float=0.0, boltzmann:bool=False, least_squares:bool=True,
                           option:str="scf", metals:tuple=spin_crossover_metals) -> np.ndarray:
    
    if interaction: interaction *= 1e3
    
    if (ms is not None) and isinstance(ms,str): ms = generate_mixed_spin_states(ls=ls, hs=hs, ms=ms)
    
    States = SpinState(states=2) if ms is None else SpinState(states=2+len(ms))
    Thermo = ThermodynamicProperties(centers=centers)

    Thermo.setref(ls=ls, hs=hs, metals=metals)

    step = extract_data(ls=ls, hs=hs, thermo=Thermo)

    States.step += [step]

    if ms is not None:
        configurations = calculate_interaction(ls=ls, hs=hs, ms=ms, interaction=interaction, states=States, thermo=Thermo, option=option)

        SlichterDrickamerModel = PhenomenologicalInteractionParameter(e_hl=step[0],
                                                                    e_ml=configurations[:,0],
                                                                    nhs=configurations[:,1],
                                                                    boltzmann_factor=configurations[:,2],
                                                                    interaction_set=configurations[:,3])
    if least_squares: boltzmann = False
    if boltzmann: least_squares = False

    if ms is not None:
        interaction = SlichterDrickamerModel.boltzmann(weights=False) if boltzmann else SlichterDrickamerModel.leastsqrs(detcoef=False)

    if temperature:
        nhs_and_gibbs  = calculate_rsm(option="direct", guess=temperature, points=points, data=States, interaction=interaction, thermo=Thermo)

        nhs_and_gibbs *= [1.0, 1e-3]

        return nhs_and_gibbs

    nhs_and_gibbs  = calculate_rsm(option="numeric", guess=298.15, points=points, data=States, interaction=interaction, thermo=Thermo)

    nhs_and_gibbs *= [1.0, 1e-3]

    return nhs_and_gibbs

#############################################################################################################################
