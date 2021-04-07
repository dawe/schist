import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

from matplotlib import cm

def get_ncolors(n_cat):
    """Return an appropriate list of colors for this dataset
    """
    
    if n_cat < 10:
        cmap = 'tab10'
    elif n_cat < 20:
        cmap = 'tab20'
    else:
        cmap = 'nipy_spectral_r'
    
    colorlist = cm.cmap_d[cmap](np.linspace(0, 1, n_cat))
    hex_list = [mpl.colors.to_hex(x) for x in colorlist]
    
    return hex_list
    
        