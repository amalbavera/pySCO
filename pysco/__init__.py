"""
The pysco package provides a user friendly data structure for thermodynamic analysis of spin crossover molecules and materials.
"""

import sys

from .read   import read
from .write  import write
from .thermo import thermo

__all__ = ["read", "write", "thermo"]

sys.tracebacklimit = 0
