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


#name='M51_V'
#kont, nx, ny = preprocess(name ,'m51.fits')
#out,scal = pethat_wavelet_scale_analysis(kont, nx, ny, 100)
#cube = wavelet_results(out,scal,name)
#cube.calc_energies()
#cube.create_gif()



#name_a = 'DD069_V'
#kont_a, nxa, nya = preprocess(name_a, 'starless_backgroundless_d69_V.fit', [[300,700],[200,600]])
#out,scal = pethat_wavelet_scale_analysis(kont_a, nxa, nya, 100)
#cube_a = wavelet_results(out,scal,name_a)
#cube_a.calc_energies()
#cube_a.create_gif()
#
#name_b = 'DD069_B'
#kont_b, nxb, nyb = preprocess(name_b, 'starless_backgroundless_d69_B.fit', [[300,700],[200,600]])
#out,scal = pethat_wavelet_scale_analysis(kont_b, nxb, nyb, 100)
#cube_b = wavelet_results(out,scal,name_b)
#cube_b.calc_energies()
#cube_b.create_gif()
#
#plot_correlation(cube_a,cube_b)




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