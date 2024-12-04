from Wavelet_Fast import *

import numpy as np
import matplotlib
from scipy import stats
from astropy.io import fits

path = str(pathlib.Path(__file__).parent.resolve())

# Testing the effect of noise in the detection of important scales

num = 300
diameter = 14
num_circles = 10
arr = np.random.standard_normal((num,num))
arr -= np.min(arr)
arr = arr / np.max(arr)
arr_rms = np.sqrt(np.var(arr))
arr_noise = arr.astype('float64')
#arr = np.zeros((num,num)) 
x, y = np.meshgrid(np.arange(num), np.arange(num))
save_array_as_FITS(arr_noise,'random_not')

num_ensemble = 20
delta = 0.1
sig_size = 15
energies = np.zeros((sig_size,50))
datacube = np.zeros((sig_size,num_ensemble,num,num))
for k in range(sig_size):
    for j in range(num_ensemble):
        datacube[k][j] = arr_noise
        x0 = np.random.randint(diameter//2,num-diameter//2,size=num_circles)
        y0 = np.random.randint(diameter//2,num-diameter//2,size=num_circles)
        for i in range(num_circles):     
            r = np.sqrt((x - x0[i])**2 + (y - y0[i])**2)
            circle = np.where(r <= diameter/ 2, (k+1)*delta*arr_rms, 0)
            datacube[k][j]+=circle
        save_array_as_FITS(datacube[k][j],'random_circles_'+str(k)+'_'+str(j))

cube1 = pethat_wavelet_scale_analysis('Random_1' ,'random_not.fits', scales_in = [2, 50, 0.8], scales_type="triplet", pixel_scale=4, distance=100000)
en_noise = cube1.calc_energies(do_plot=False)
for i in range(sig_size):
    for j in range(num_ensemble):
        cube = pethat_wavelet_scale_analysis('Random_Circle_'+str(i)+'_'+str(j) ,'random_circles_'+str(i)+'_'+str(j)+'.fits', scales_in = [2, 50, 0.8], scales_type="triplet", pixel_scale=4, distance=100000)
        energies[i] += cube.calc_energies(do_plot=False) #- en_noise
energies = energies/num_ensemble

def update(frame):
    img.set_data(cube1.scales,energies[frame])
    tx.set_text('Energy of the data with circles minus the energy of noise, scale = '+str(np.round((frame+1)*delta,1))+" sigma")

fig , ax = plt.subplots(1,1,dpi=300)
img, = ax.plot(cube1.scales,energies[0], animated=True)
ax.set_xlabel(r'Scale (pixels)', fontsize=11)
ax.set_ylabel(r'Wavelet Energy of random circles', fontsize=11)
ax.set_ylim([np.min(energies),np.max(energies)+1])
tx = ax.set_title('Energy of circles + noise, scale = '+str(np.round(delta,1))+" sigma")
ani = animation.FuncAnimation(fig=fig, func=update, frames=sig_size, interval=200)
ani.save(filename=path+"\\Output\\Random_circle_noise_energies_animation_20.gif", writer='ffmpeg',codec="libx264")
