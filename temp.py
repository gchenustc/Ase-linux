import numpy as np
import ase
from ase.io import read
from ase.io import write
from ase.spacegroup import crystal
from ase.calculators.vasp import Vasp
import matplotlib.pyplot as plt
import matplotlib
from ase.dft.dos import DOS


import matplotlib
print(matplotlib.matplotlib_fname())