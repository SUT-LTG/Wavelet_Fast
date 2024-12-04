from Wavelet_Fast import *

import numpy as np
import matplotlib
from scipy import stats
from astropy.io import fits

path = str(pathlib.Path(__file__).parent.resolve())

cube = pethat_wavelet_scale_analysis('M51_V' ,'m51.fits', scales_in = [2, 100, 1], scales_type="triplet", pixel_scale=4, distance=100000)
cube.save_layers(unit='pixels')
cube.save_FITS()
cube.calc_energies(unit='pixels')
cube.create_gif(unit='pixels')
