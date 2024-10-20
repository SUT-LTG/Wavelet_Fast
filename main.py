from Wavelet_Fast import *

import numpy as np
import matplotlib
from scipy import  fft
from astropy.io import fits

#------------------------------
#--------- test data ----------
#------------------------------

#kont, outname, nx,ny = preprocess('DDO168_V','starless_backgroundless_d168.fit')

#kont, outname, nx,ny = preprocess('M51_V','m51.fits')

kont, outname, nx,ny = preprocess('DD069_V','starless_backgroundless_d69_V.fit',[[300,700],[200,600]])

#kont, outname, nx,ny = preprocess('DD069_B','starless_backgroundless_d69_B.fit',[[300,700],[200,600]])


energies, en_scales, outputs = pethat_wavelet_scale_analysis(kont,outname,nx,ny, 100)

#--------------------------------------------------------------------------------------------------------------------------------------------------

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