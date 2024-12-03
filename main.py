from Wavelet_Fast import *

import numpy as np
import matplotlib
from scipy import stats
from astropy.io import fits

path = str(pathlib.Path(__file__).parent.resolve())
#------------------------------
#--------- test data ----------
#------------------------------

#kont, outname, nx, ny = preprocess('DDO168_V', 'starless_backgroundless_d168.fit')

"""First Test"""
#cube = pethat_wavelet_scale_analysis('M51_V' ,'m51.fits', scales_in = [2, 100, 1], scales_type="triplet", pixel_scale=4, distance=100000)
##cube.save_layers(unit='pixels')
#cube.save_FITS()
#cube.calc_energies(unit='pixels')
#cube.create_gif(unit='pixels')


"""Testing the cross-correlation"""
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

#cube = pethat_wavelet_scale_analysis('M74_Unknown', 'starless_m74_unknown.fit', scales_in = [5, 50, 10], scales_type="triplet")
#cube.save_layers()
#cube.calc_energies()
#cube.create_gif()


#(nx,ny) = (300,300)
#name = 'Random'
#outname = 'Output/'+name+'/'+name+'_'
#kont = np.random.rand(nx,ny)

#plt.imshow(kont)
#plt.colorbar()
#plt.savefig(path+"\\"+outname +'original.png', dpi=200)
##plt.show()

"""--------------------------------------------------------------------------------------------------------------------------------------------------"""
# Testing the effect of noise in the detection of important scales

"""
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

"""

# Same test with Gaussian stars

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
delta = 0.5
sig_size = 20
energies = np.zeros((sig_size,70))
datacube = np.zeros((sig_size,num_ensemble,num,num))
for k in range(sig_size):
    for j in range(num_ensemble):
        datacube[k][j] = arr_noise
        x0 = np.random.randint(diameter//2,num-diameter//2,size=num_circles)
        y0 = np.random.randint(diameter//2,num-diameter//2,size=num_circles)
        for i in range(num_circles):     
            r = np.sqrt((x - x0[i])**2 + (y - y0[i])**2)
            circle = (k+1)*delta*arr_rms*np.exp(-r**2*2/diameter**2)
            datacube[k][j]+=circle
        save_array_as_FITS(datacube[k][j],'random_circles_'+str(k)+'_'+str(j))

cube1 = pethat_wavelet_scale_analysis('Random_1' ,'random_not.fits', scales_in = [5, 70, 0.1], scales_type="triplet")
en_noise = cube1.calc_energies(do_plot=False)
for i in range(sig_size):
    for j in range(num_ensemble):
        cube = pethat_wavelet_scale_analysis('Random_Circle_'+str(i)+'_'+str(j) ,'random_circles_'+str(i)+'_'+str(j)+'.fits', scales_in = [3, 70, 0.25], scales_type="triplet")
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
deriv = (energies[:,1:69]-energies[:,0:68])/0.25
def update2(frame):
    img.set_data(cube1.scales[1:69],deriv[frame])
    tx.set_text('Derivative of Energy of the data with stars minus the energy of noise, scale = '+str(np.round((frame+1)*delta,1))+" sigma")

fig , ax = plt.subplots(1,1,dpi=300)
img, = ax.plot(cube1.scales[1:69],deriv[0], animated=True)
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
    for i in range(67):
        if deriv[j,i]>0 and deriv[j,i+1]<0 and mem==0:
            peak[j] = 5+i*0.1
            mem=1

plt.clf()
plt.scatter(peak,snr)
plt.xlabel(r'location of the Peak (pixels)', fontsize=11)
plt.ylabel(r'SNR', fontsize=11)
plt.savefig(path+"\\Output\\"+'peak_to_snr_D6.png', dpi=200)

pear = stats.pearsonr(peak,snr)
print('Pearson correlation coefficient is: '+str(pear.statistic)+'. The p-value is: '+str(pear.pvalue) )


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