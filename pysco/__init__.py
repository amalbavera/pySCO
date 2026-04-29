"""
pySCO provides a user friendly data structure for
thermodynamic analysis of spin crossover molecules and materials.

Modules
-------
    spin_crossover_energy
        Compute the spin-crossover energy between the high- and low-spin states.
    
    transition_temperature
        Compute the transition temperature.

    high_spin_population
        Compute the thermal variation of the relative high-spin population, or magnetization.

    gibbs_free_energy
        Compute the isothermal variation of the Gibbs free energy.

Output Parsers
--------------
    read_vasp
        Read output files generated with the Vienna Ab initio Simulation Program.

    read_espresso
        Read output files generated with Quantum Espresso.

    read_gaussian
        Read output files generated with Gaussian.

    read_orca
        Read output files generated with Orca.

    read_nwchem
        Read output files generated with NWChem.

    read_pyscf
        Read output objects and dictionaries generated with pySCF.

Thermodynamic Models
--------------------
    RegularSolution
        Use the simple regular solution to model the spin conversion.

    SlichterDrickamer
        Use the extended Slichter and Drickamer model to model the spin conversion by including a
        phenomenological interaction parameter.

Microscopic Models
------------------
    NoneSoFar
        No microscopic models are available at the time.
"""
#
### All modules available
#
__all__ = ["read_vasp",
           "read_gaussian", "read_nwchem", "read_orca", "read_pyscf",
           "models", "read",
           "RegularSolution", "SlichterDrickamer"]
#
### Base-line objects
#
from .io      import read
from .models  import models
#
### Materials science codes
#
from .io.read import Vasp     as read_vasp
#
### Quantum chemistry codes
#
from .io.read import Gaussian as read_gaussian
from .io.read import NWChem   as read_nwchem
from .io.read import Orca     as read_orca
from .io.read import pySCF    as read_pyscf
#
### Thermodynamic models
#
from .models.thermodynamic import RegularSolution
from .models.thermodynamic import SlichterDrickamer
#
### Microscopic models
#
#from .models.microscopic   import 

