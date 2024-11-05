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


#cube = pethat_wavelet_scale_analysis('M51_V' ,'m51.fits', start = 2, step_length = 1, n = 100)
#cube.save_layers()
#cube.save_FITS()
#cube.calc_energies()
#cube.create_gif()



cube_a = pethat_wavelet_scale_analysis('DD069_V', 'starless_backgroundless_d69_V.fit', [[300,700],[200,600]], 2, 1, 100)
cube_a.save_layers()
cube_a.save_FITS()
cube_a.calc_energies()
cube_a.create_gif()

cube_b = pethat_wavelet_scale_analysis('DD069_B', 'starless_backgroundless_d69_B.fit', [[300,700],[200,600]], 2, 1, 100)
cube_b.save_layers()
cube_b.save_FITS()
cube_b.calc_energies()
cube_b.create_gif()

plot_correlation(cube_a,cube_b)




#(nx,ny) = (300,300)
#name = 'Random'
#outname = 'Output/'+name+'/'+name+'_'
#kont = np.random.rand(nx,ny)

#plt.imshow(kont)
#plt.colorbar()
#plt.savefig(path+"\\"+outname +'original.png', dpi=200)
##plt.show()









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