"""
The pySCO package provides a user friendly data structure for
thermodynamic analysis of spin crossover molecules and materials.
"""

import sys

from .read   import read
from .thermo import thermo

__all__ = ["read", "thermo"]

sys.tracebacklimit = 0

