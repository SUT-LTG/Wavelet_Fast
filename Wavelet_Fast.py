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
#---------------------- Pre-Processing ---------------------------
#-----------------------------------------------------------------


def preprocess(name,filename,crop_cor = []):
    path = str(pathlib.Path(__file__).parent.resolve())
    #------------------------------
    #-- Retrieving observational data
    #------------------------------
    outname = 'Output/'+name+'/'+name+'_'
    adata = fits.open(path+'\\data\\'+filename)
    a = adata[0].data
    #------------------------------
    #-- Cropping and zero padding
    #------------------------------
    if crop_cor != []:
        a=a[crop_cor[0][0]:crop_cor[0][1],crop_cor[1][0]:crop_cor[1][1]]

    kont = np.array(a.astype(float))
    kont = kont/np.max(kont)  # Normalize to 1

    plt.imshow(kont)
    plt.colorbar()
    plt.savefig(path+"\\"+outname +'original.png', dpi=200)
    plt.show()

    ny = kont[0,:].size
    nx = kont[:,0].size
    kont = np.pad(kont,[(0,int(np.max([-nx+ny,0]))), (0,int(np.max([nx-ny,0])))], mode='constant') 

    return kont, outname, nx, ny # Returns the final map and the save path

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

#-----------------------------------------------------------------
#--------------------- Wavelet Scale Analysis --------------------
#-----------------------------------------------------------------

def pethat_wavelet_scale_analysis(kont, outname, nx, ny, n = 100, dosum = False):

    path = str(pathlib.Path(__file__).parent.resolve())

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
            if dosum == True:
                sums += map1 
            plt.imshow(map1, origin='lower', interpolation='nearest')
            plt.colorbar()
            plt.title('PetHat wavelet, scale='+str(np.round(en_scales[i],5)), fontsize=16)
            plt.savefig(path+"\\"+outname +'wavelet_pethat_fast_'+str(i+1)+'.png', dpi=200)
            plt.clf()
            bar()

    if dosum == True:
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


    return energies, en_scales, wavelet_coeffs


#-----------------------------------------------------------------
#-------------------- Wavelet Scale Correlation ------------------
#-----------------------------------------------------------------


def pethat_wavelet_scale_correlation(kont_a, outname_a, kont_b, outname_b, nx, ny, n = 100):    
    
    path = str(pathlib.Path(__file__).parent.resolve())
    energies_a, en_scales_a, outputs_a = pethat_wavelet_scale_analysis(kont_a, outname_a, nx, ny, n)
    energies_b, en_scales_b, outputs_b = pethat_wavelet_scale_analysis(kont_b, outname_b, nx, ny, n)

    print("The normal correlation between the maps is: ", (np.sum( kont_a * np.conjugate(kont_b) )/np.sqrt( np.sum(np.abs(kont_a)**2) * np.sum(np.abs(kont_b)**2) )))
    
    corr = np.zeros(n) #Calculating the cross-correlation
    for i in range(n):
        corr[i] = np.real(np.sum(outputs_a[i]*np.conjugate(outputs_b[i]))/np.sqrt(energies_a[i]*energies_b[i]))
    corr_err = np.sqrt(1 - corr**2)/ np.sqrt(nx*ny/en_scales_a**2 - 2)

    plt.errorbar(en_scales_a,corr,yerr = corr_err,fmt = '.')
    plt.xlabel(r'Scale', fontsize=14)
    plt.ylabel(r'Correlation of B and V', fontsize=14)
    plt.title('Correlation to scale', fontsize=16)
    plt.savefig(path+'\\'+'Output/'+'6cm_pethat_corr_smooth.png', dpi=400)
    plt.show()
    
    return
      

class data_cube:
    def __init__(self, cube, scales , name):
        self.cube = cube
        self.scales = scales
        self.name = name
        self.outname = 'Output/'+name+'/'+name+'_'  

    def create_gif(self):
        def update(frame,cube,scales,name):
            img.set_array(cube[frame])
            tx.set_text('PetHat wavelet animation of ' + name + ' Scale = '+ str(np.round(scales[frame],5)))
        
        path = str(pathlib.Path(__file__).parent.resolve())
        
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(1,1,1)
        img = ax.imshow(self.cube[0], origin='lower', interpolation='nearest',vmin=np.min(self.cube), vmax=np.max(self.cube), animated=True)
        cb = fig.colorbar(img)
        tx = ax.set_title('PetHat wavelet animation of ' + self.name + ' Scale = '+ str(np.round(self.scales[0],5)), fontsize=16)

        ani = animation.FuncAnimation(fig=fig, func=update, fargs=(self.cube,self.scales,self.name), frames=len(self.cube)-1, interval=100)
        ani.save(filename=path+"\\"+self.outname +'animation.gif', writer="ffmpeg",codec="libx264")
        
        return
    
    def save_FITS(self):
        path = str(pathlib.Path(__file__).parent.resolve())
        hdu_new = fits.PrimaryHDU(np.real(self.cube))
        hdu_new.writeto(path+'\\'+self.outname+'wavelet_coeffs.fits',overwrite=True)
        return
        