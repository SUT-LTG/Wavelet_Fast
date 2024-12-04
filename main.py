from Wavelet_Fast import *

import numpy as np
import matplotlib
from scipy import stats
from astropy.io import fits

path = str(pathlib.Path(__file__).parent.resolve())
#------------------------------
#--------- test data ----------
#------------------------------


"""Testing large data"""
#cube = pethat_wavelet_scale_analysis('M74_IR', 'JWST-M74-f444w_fixed.fits', scales_in = [10, 50, 8], scales_type="triplet")
#cube.save_layers()
#cube.save_FITS()
#cube.calc_energies()
#cube.create_gif()

"""Testing data of M74"""
#cube = pethat_wavelet_scale_analysis('M74_Amateur', 'starless_m74_ama.fit', scales_in = [5, 150, 1], scales_type="triplet")
#cube.save_layers()
#cube.calc_energies()
#cube.create_gif()

#cube = pethat_wavelet_scale_analysis('M74_g', 'starless_m74_unknown.fit', scales_in = [5, 50, 10], scales_type="triplet")
#cube.save_layers()
#cube.calc_energies()
#cube.create_gif()


#--------------------------------------------------------------------------------------------------------------------------------------------------
" Some old tests " 

#nx = 1000
#ny = 1000
#coeff = 2*np.pi/(2*nx-1)
#sca = 50 * coeff

#pet_r = scipy.fft.fftshift(scipy.fft.fft2(scipy.fft.ifftshift(pethat_phi_func(np.zeros((2*nx-1,2*ny-1)),sca))))
#pet_r_rad = np.real(pet_r)
#pet_r_rad = pet_r_rad[int(np.max([nx,ny]))]

#plt.imshow(np.abs(np.real(pet_r)))
#plt.colorbar()
#plt.show()

#plt.plot(np.arange(-int(np.max([nx,ny])),int(np.max([nx,ny]))-1),pet_r_rad)
#plt.show()