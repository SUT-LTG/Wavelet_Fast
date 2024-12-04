from Wavelet_Fast import *

import numpy as np
import matplotlib
from scipy import stats
from astropy.io import fits

path = str(pathlib.Path(__file__).parent.resolve())

num = 300
diameter = 6
num_circles = 25
arr = np.random.standard_normal((num,num))
arr -= np.min(arr)
arr = arr / np.max(arr)
arr_rms = np.sqrt(np.var(arr))
arr_noise = arr.astype('float64')
#arr = np.zeros((num,num)) 
x, y = np.meshgrid(np.arange(num), np.arange(num))
save_array_as_FITS(arr_noise,'random_not')

num_ensemble = 20
delta = 0.3
sig_size = 25
energies = np.zeros((sig_size,50))
datacube = np.zeros((sig_size,num_ensemble,num,num))
for k in range(sig_size):
    for j in range(num_ensemble):
        datacube[k][j] = arr_noise
        x0 = np.random.randint(diameter*2,num-diameter*2,size=num_circles)
        y0 = np.random.randint(diameter*2,num-diameter*2,size=num_circles)
        for i in range(num_circles):     
            r = np.sqrt((x - x0[i])**2 + (y - y0[i])**2)
            circle = (k+1)*delta*arr_rms*np.exp(-r**2*2/diameter**2)
            datacube[k][j]+=circle
        save_array_as_FITS(datacube[k][j],'random_circles_'+str(k)+'_'+str(j))

scales_triplet = [7, 50, 0.1]
cube1 = pethat_wavelet_scale_analysis('Random_1' ,'random_not.fits', scales_in = scales_triplet, scales_type="triplet")
en_noise = cube1.calc_energies(do_plot=False)
for i in range(sig_size):
    for j in range(num_ensemble):
        cube = pethat_wavelet_scale_analysis('Random_Circle_'+str(i)+'_'+str(j) ,'random_circles_'+str(i)+'_'+str(j)+'.fits', scales_in = scales_triplet, scales_type="triplet")
        energies[i] += cube.calc_energies(do_plot=False) - en_noise
energies = energies/num_ensemble

def update(frame):
    img.set_data(cube1.scales,energies[frame])
    tx.set_text('Energy of the data with stars minus the energy of noise, scale = '+str(np.round((frame+1)*delta,1))+" sigma")

fig , ax = plt.subplots(1,1,dpi=300)
img, = ax.plot(cube1.scales,energies[0], animated=True)
ax.set_xlabel(r'Scale (pixels)', fontsize=11)
ax.set_ylabel(r'Wavelet Energy of random circles', fontsize=11)
ax.set_ylim([np.min(energies),np.max(energies)+1])
tx = ax.set_title('Energy of the data with stars minus the energy of noise, scale = '+str(np.round(delta,1))+" sigma")
ani = animation.FuncAnimation(fig=fig, func=update, frames=sig_size, interval=200)
ani.save(filename=path+"\\Output\\Random_stars_energies_minus_noise_animation_D6.gif", writer='ffmpeg',codec="libx264")

plt.clf()
deriv = (energies[:,1:49]-energies[:,0:48])/0.1
def update2(frame):
    img.set_data(cube1.scales[1:49],deriv[frame])
    tx.set_text('Derivative of Energy of the data with stars minus the energy of noise, scale = '+str(np.round((frame+1)*delta,1))+" sigma")

fig , ax = plt.subplots(1,1,dpi=300)
img, = ax.plot(cube1.scales[1:49],deriv[0], animated=True)
plt.xlabel(r'Scale (pixels)', fontsize=11)
plt.ylabel(r'Derivative of Energy', fontsize=11)
ax.set_ylim([np.min(deriv)-1,np.max(deriv)+1])
tx = ax.set_title('Derivative of Energy of the data with stars minus the energy of noise, scale = '+str(np.round(delta,1))+" sigma")
ani = animation.FuncAnimation(fig=fig, func=update2, frames=sig_size, interval=200)
ani.save(filename=path+"\\Output\\Random_stars_energies_derivative_animation_D6.gif", writer='ffmpeg',codec="libx264")

peak = np.zeros(sig_size)
snr = np.zeros(sig_size)
for j in range(sig_size):
    snr[j] = (j+1)*delta
    mem = 0
    for i in range(47):
        if deriv[j,i]>0 and deriv[j,i+1]<0 and mem==0:
            peak[j] = 7+(i+1)*0.1
            mem = 1
    if mem == 0:
        peak[j] = np.nan

plt.clf()
plt.scatter(peak,snr)
plt.title(r'location of the Peak vs the SNR')
plt.xlabel(r'location of the Peak (pixels)', fontsize=11)
plt.ylabel(r'SNR', fontsize=11)
plt.savefig(path+"\\Output\\"+'peak_to_snr_D6.png', dpi=200)

spear = stats.spearmanr(peak,snr)
print('Spearman correlation coefficient and the associated p-value are: '+str(spear))

pear = stats.pearsonr(peak,snr)
print('Pearson correlation coefficient and the associated p-value are: '+str(pear))
