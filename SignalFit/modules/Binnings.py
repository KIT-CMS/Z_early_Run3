import numpy as np
from collections import OrderedDict

"""
mass bins for rebinning mT and fits
"""

mass_bins_w = np.linspace(0, 120, num=20+1)
mass_bins_w_ptOverMt = np.linspace(0, 3, num=30+1)
mass_bins_z = np.linspace(60, 120, num=30+1)

mass_bins_w_scan = OrderedDict()
for i in range(1):
    mass_bins_w_scan[i] = mass_bins_w[i:]
