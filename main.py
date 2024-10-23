from Wavelet_Fast import *

import numpy as np
import matplotlib
from scipy import  fft
from astropy.io import fits

path = str(pathlib.Path(__file__).parent.resolve())
#------------------------------
#--------- test data ----------
#------------------------------

#kont, outname, nx, ny = preprocess('DDO168_V', 'starless_backgroundless_d168.fit')

#kont, outname, nx, ny = preprocess('M51_V' ,'m51.fits')

#kont_a, outname_a, nx, ny = preprocess('DD069_V', 'starless_backgroundless_d69_V.fit', [[300,700],[200,600]])

#kont_b, outname_b, nx, ny = preprocess('DD069_B', 'starless_backgroundless_d69_B.fit', [[300,700],[200,600]])

(nx,ny) = (300,300)
name = 'Random'
outname = 'Output/'+name+'/'+name+'_'
kont = np.random.rand(nx,ny)

plt.imshow(kont)
plt.colorbar()
plt.savefig(path+"\\"+outname +'original.png', dpi=200)
#plt.show()

energies, en_scales, outputs = pethat_wavelet_scale_analysis(kont, outname, nx, ny, 100)

cube = data_cube(np.real(outputs),en_scales,name)
cube.create_gif()
cube.save_FITS()

#pethat_wavelet_scale_correlation(kont_a, outname_a, kont_b, outname_b, nx, ny, 100)

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