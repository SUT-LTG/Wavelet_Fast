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


#cube = pethat_wavelet_scale_analysis('M51_V' ,'m51.fits', scales_in = [2, 100, 1], scales_type="triplet", pixel_scale=4, distance=100000)
##cube.save_layers(unit='pixels')
#cube.save_FITS()
#cube.calc_energies(unit='pixels')
#cube.create_gif(unit='pixels')
#
#
#
#cube_a = pethat_wavelet_scale_analysis('DD069_V', 'starless_backgroundless_d69_V.fit', [[300,700],[200,600]], scales_in = [2, 100, 1], scales_type="triplet", pixel_scale=1.134, distance=800000)
##cube_a.save_layers()
#cube_a.save_FITS()
#cube_a.calc_energies()
#cube_a.create_gif()
#
#cube_b = pethat_wavelet_scale_analysis('DD069_B', 'starless_backgroundless_d69_B.fit', [[300,700],[200,600]],  scales_in = [2, 100, 1],scales_type="triplet", pixel_scale=1.134, distance=800000)
##cube_b.save_layers()
#cube_b.save_FITS()
#cube_b.calc_energies()
#cube_b.create_gif()
##
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
# Testing some properties 
"""
num = 300
lamb = 10
diameter = 14
num_circles = 10
det_lim = 0.2
arr = np.random.poisson(lamb,(num,num))
arr_rms = np.sqrt(np.var(arr))
arr = arr.astype('float64')
#arr = np.zeros((num,num)) 
x, y = np.meshgrid(np.arange(num), np.arange(num))

save_array_as_FITS(arr,'random_not')
x0 = np.random.randint(diameter,num-diameter,size=num_circles)
y0 = np.random.randint(diameter,num-diameter,size=num_circles)
for i in range(num_circles):
    r = np.sqrt((x - x0[i])**2 + (y - y0[i])**2)
    circle = np.where(r <= diameter/ 2, det_lim*arr_rms, 0)
    #circle = det_lim*arr_rms*np.exp(-2*(r/diameter)**2)
    arr+=circle

save_array_as_FITS(arr,'random_circles')

cube = pethat_wavelet_scale_analysis('Random_Circle' ,'random_circles.fits', scales_in = [2, 80, 0.5], scales_type="triplet", pixel_scale=4, distance=100000)
cube1 = pethat_wavelet_scale_analysis('Random_1' ,'random_not.fits', scales_in = [2, 80, 0.5], scales_type="triplet", pixel_scale=4, distance=100000)
#cube.save_layers(unit='pixels')
#cube.save_FITS()
en = cube.calc_energies(unit='pixels')
en1 = cube1.calc_energies(unit='pixels',do_plot=False)
#cube.create_gif(unit='pixels')
plt.plot(cube.scales,en-en1)
plt.xlabel(r'Scale (pixels)', fontsize=11)
plt.ylabel(r'Wavelet Energy Difference ', fontsize=11)
plt.title('PetHat Wavelet Energies - noise + circles R='+str(diameter)+' - only noise', fontsize=11)
plt.savefig(path+'\\Output\\pethat_energy_smooth_circles_'+str(det_lim)+'sig.png', dpi=300)
#plt.show()

"""

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