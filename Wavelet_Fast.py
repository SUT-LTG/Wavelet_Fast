import numpy as np
import scipy, time
import matplotlib
import pathlib

from alive_progress import alive_bar
import matplotlib.animation as animation

from scipy import integrate, fft ,signal
import pylab as plt
from scipy import io
from astropy.io import fits
from scipy import stats

#-----------------------------------------------------------------
#---------------------- The Pet-Hat function ---------------------
#-----------------------------------------------------------------

def pethat_k(k):
    if ((np.pi <= k)and(k <= (4.0*np.pi))): 
        return np.cos(np.pi * 0.5 * np.log2(k * 0.5/np.pi))**2
    else:
        return 0.


def pethat_phi_func(data,a_scale):
    ny = data[0,:].size
    nx = data[:,0].size
    val_map = np.zeros((nx, ny), np.float32)
    for i in np.arange(0, nx):
        for j in np.arange(0, ny):
            val_map[i][j] = pethat_k(np.sqrt((i-nx/2)**2 + (j-ny/2)**2)*a_scale)
    return val_map

def pethat_wavelet_scale_analysis(kont,outname, n = 100):

    fft_shape = 2*np.array(kont.shape) - 1
    fft_kont = scipy.fft.fftshift(scipy.fft.fft2(kont,fft_shape))
    fnx = fft_kont[:,0].size
    coeff = 2*np.pi/fnx 

    energies = np.zeros(n)
    en_scales = 2 + np.arange(0,n) / n * 100
    sums = np.zeros((nx,ny))
    wavelet_coeffs = np.clongdouble(np.zeros((n,nx,ny)))

    with alive_bar(n) as bar:
        for i in range(n):
            a = en_scales[i] * coeff
            map1 = scipy.fft.ifft2(scipy.fft.ifftshift(fft_kont*pethat_phi_func(fft_kont,a)))
            map1 = map1[0:nx,0:ny]
            energies[i] = np.sum(np.abs(map1)**2)
            wavelet_coeffs[i] = map1
            map1 = np.real(map1)
            #map1 = map1 * np.heaviside(map1,0)
            if (i<=80 and 2<i) or True:
                sums += map1 
            plt.imshow(map1, origin='lower', interpolation='nearest')
            plt.colorbar()
            plt.title('PetHat wavelet, scale='+str(np.round(en_scales[i],5)), fontsize=16)
            plt.savefig(path+"\\"+outname +'wavelet_pethat_fast_'+str(i+1)+'.png', dpi=200)
            plt.clf()
            bar()
    return sums ,energies, en_scales, wavelet_coeffs


path = str(pathlib.Path(__file__).parent.resolve())
#------------------------------
#-- observational data
#------------------------------

outname = 'Output/DDO168_V/DDO168_V_'
adata = fits.open(path+'\\data\\'+'starless_backgroundless_d168.fit')

#outname = 'Output/M51_V/M51_V_'
#adata = fits.open(path+'\\data\\'+'m51.fits')

#outname = 'Output/DD069_V/DD069_V_'
#adata = fits.open(path+'\\data\\'+'starless_backgroundless_d69_V.fit')

a = adata[0].data

#------------------------------
#-- cropping and zero padding
#------------------------------

#a=a[300:700,200:600]

kont = np.array(a.astype(float))
kont = kont/np.max(kont)

plt.imshow(kont)
plt.colorbar()
plt.show()

ny = kont[0,:].size
nx = kont[:,0].size

print(kont.shape,np.max(kont),np.min(kont))
kont = np.pad(kont,[(0,int(np.max([-nx+ny,0]))), (0,int(np.max([nx-ny,0])))], mode='constant')
print(kont.shape,np.max(kont),np.min(kont))

#--------------------------------------------------------------------------------------------------------------------------------------------------
"""
outname_B = 'Output/DD069_B/DD069_B_'
adata_B = fits.open(path+'\\data\\'+'starless_backgroundless_d69_B.fit')

b = adata_B[0].data


b=b[300:700,200:600]

kont_b = np.array(b.astype(float))
kont_b = kont_b/np.max(kont_b)


plt.imshow(kont_b)
plt.colorbar()
plt.show()

nyb = kont_b[0,:].size
nxb = kont_b[:,0].size

print(kont_b.shape,np.max(kont_b),np.min(kont_b))
kont_b = np.pad(kont_b,[(0,int(np.max([-nxb+nyb,0]))), (0,int(np.max([nxb-nyb,0])))], mode='constant')
print(kont_b.shape,np.max(kont_b),np.min(kont_b))
"""
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

sums ,energies, en_scales, outputs = pethat_wavelet_scale_analysis(kont,outname, 100)

"""sums_b ,energies_b, en_scales_b, outputs_b = pethat_wavelet_scale_analysis(kont_b,outname_B, 100)

n=en_scales_b.size
print(np.sum(kont*np.conjugate(kont_b))/np.sqrt(np.sum(np.abs(kont)**2)*np.sum(np.abs(np.conjugate(kont_b))**2)))
corr = np.zeros(n)
for i in range(n):
    corr[i] = np.real(np.sum(outputs[i]*np.conjugate(outputs_b[i]))/np.sqrt(energies[i]*energies_b[i]))


plt.plot(en_scales,corr)
plt.xlabel(r'Scale', fontsize=14)
plt.ylabel(r'Correlation of B and V', fontsize=14)
plt.title('Correlation to scale', fontsize=16)
plt.savefig(path+'\\'+'Output/'+'6cm_pethat_corr_smooth.png', dpi=200)
plt.show()
"""

plt.imshow(sums, origin='lower', interpolation='nearest')
plt.colorbar()
plt.title('Sum of the PetHat Wavelet Coefficients', fontsize=16)
plt.savefig(path+"\\"+outname +'wavelet_pethat_fast_total.png', dpi=200)
#plt.show()
plt.clf()


plt.plot(en_scales,energies)
plt.xlabel(r'Scale', fontsize=14)
plt.ylabel(r'Wavelet Energy', fontsize=14)
plt.title('PetHat Wavelet Energies', fontsize=16)
plt.savefig(path+'\\'+outname+'6cm_pethat_energy_smooth.png', dpi=200)
#plt.show()
plt.clf()

"""
plt.imshow(sums_b, origin='lower', interpolation='nearest')
plt.colorbar()
plt.title('Sum of the PetHat Wavelet Coefficients', fontsize=16)
plt.savefig(path+"\\"+outname_B +'wavelet_pethat_fast_total.png', dpi=200)
#plt.show()
plt.clf()

plt.plot(en_scales_b,energies_b)
plt.xlabel(r'Scale', fontsize=14)
plt.ylabel(r'Wavelet Energy', fontsize=14)
plt.title('PetHat Wavelet Energies', fontsize=16)
plt.savefig(path+'\\'+outname_B+'6cm_pethat_energy_smooth.png', dpi=200)
#plt.show()
"""
hdu_new = fits.PrimaryHDU(np.real(outputs))
hdu_new.writeto(path+'\\'+outname+'wavelet_coeffs.fits',overwrite=True)