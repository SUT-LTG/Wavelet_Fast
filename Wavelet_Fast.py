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
path = str(pathlib.Path(__file__).parent.resolve())

def preprocess(name,filename,crop_cor = []):
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
    #plt.show()
    plt.clf()

    ny = kont[0,:].size
    nx = kont[:,0].size
    kont = np.pad(kont,[(0,int(np.max([-nx+ny,0]))), (0,int(np.max([nx-ny,0])))], mode='constant') 

    return kont, nx, ny # Returns the final map and the save path

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


def pethat_wavelet_scale_analysis(kont, nx, ny, n = 100):

    fft_shape = 2*np.array(kont.shape) - 1
    fft_kont = scipy.fft.fftshift(scipy.fft.fft2(kont,fft_shape))
    fnx = fft_kont[:,0].size
    coeff = 2*np.pi/fnx 

    en_scales = 2 + np.arange(0,n) / n * 100
    wavelet_coeffs = np.clongdouble(np.zeros((n,nx,ny)))

    with alive_bar(n) as bar: 
        for i in range(n):
            a = en_scales[i] * coeff
            map1 = scipy.fft.ifft2(scipy.fft.ifftshift(fft_kont*pethat_phi_func(fft_kont,a)))
            map1 = map1[0:nx,0:ny]
            wavelet_coeffs[i] = map1
            bar()

    return wavelet_coeffs,en_scales

#-----------------------------------------------------------------
#------------------ Wavelet Results as an Object -----------------
#-----------------------------------------------------------------

      
class wavelet_results:
    def __init__(self, cube, scales , name):
        self.cube = cube
        self.scales = scales
        self.name = name
        self.outname = 'Output/'+name+'/'+name+'_'  

    def create_gif(self):
        def update(frame,cube,scales,name):
            img.set_array(cube[frame])
            tx.set_text('PetHat wavelet animation of ' + name + ' Scale = '+ str(np.round(scales[frame],5)))
        
        our_data = np.real(self.cube)
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(1,1,1)
        img = ax.imshow(our_data[0], origin='lower', interpolation='nearest',vmin=np.min(our_data), vmax=np.max(our_data), animated=True)
        cb = fig.colorbar(img)
        tx = ax.set_title('PetHat wavelet animation of ' + self.name + ' Scale = '+ str(np.round(self.scales[0],5)), fontsize=16)

        ani = animation.FuncAnimation(fig=fig, func=update, fargs=(our_data,self.scales,self.name), frames=len(self.cube)-1, interval=100)
        ani.save(filename=path+"\\"+self.outname +'animation.gif', writer="ffmpeg",codec="libx264")
        
        return
    
    def save_FITS(self):
        hdu_new = fits.PrimaryHDU(np.real(self.cube))
        hdu_new.writeto(path+'\\'+self.outname+'wavelet_coeffs.fits',overwrite=True)
        return
        
    def save_layers(self):
        for i in range(len(self.cube)):
            plt.imshow(np.real(self.cube[i]), origin='lower', interpolation='nearest')
            plt.colorbar()
            plt.title('PetHat wavelet coefficient of '+self.name+', scale='+str(np.round(self.scales[i],5))+' pixels', fontsize=16)
            plt.savefig(path+"\\"+self.outname +'wavelet_pethat_fast_'+str(i+1)+'.png', dpi=300)
            plt.clf()
        return

    def calc_energies(self, do_plot=True):
        energies = np.zeros(len(self.cube))
        for i in range(len(self.cube)):
            energies[i] = np.sum(np.abs(self.cube[i])**2)
        if do_plot:
            plt.clf()
            plt.plot(self.scales,energies)
            plt.xlabel(r'Scale (pixels)', fontsize=14)
            plt.ylabel(r'Wavelet Energy of '+self.name, fontsize=14)
            plt.title('PetHat Wavelet Energies', fontsize=16)
            plt.savefig(path+'\\'+self.outname+'pethat_energy_smooth.png', dpi=300)
            #plt.show()
            plt.clf()
        return energies
    
    def calc_sum(self, do_plot=True):
        sums = np.sum(np.real(self.cube),axis=0)
        if do_plot:
            plt.imshow(sums, origin='lower', interpolation='nearest')
            plt.colorbar()
            plt.title('Sum of the PetHat Wavelet Coefficients', fontsize=16)
            plt.savefig(path+"\\"+self.outname +'wavelet_pethat_fast_total.png', dpi=200)
            #plt.show()
            plt.clf()
        return sums


#-----------------------------------------------------------------
#-------------------- Wavelet Scale Correlation ------------------
#-----------------------------------------------------------------


def give_names(str_a,str_b):
    name = ''
    for i in range(len(str_a)):
        if str_a[i] == '_':
            str_a = str_a[i+1:]
            str_b = str_b[i+1:]
            break
        name+=(str_a[i])  
    return name , str_a, str_b

def plot_correlation(cube_a, cube_b):
    energies_a = cube_a.calc_energies(do_plot=False)
    energies_b = cube_b.calc_energies(do_plot=False)
    if (not np.allclose(cube_a.scales,cube_b.scales)) or np.shape(cube_a.cube) != np.shape(cube_b.cube):
        print('Error: Cubes are not compatible.')
        return
    (n,nx,ny) = np.shape(cube_a.cube)
    corr = np.zeros(n) #Calculating the cross-correlation and its error
    for i in range(n):
        corr[i] = np.real(np.sum(cube_a.cube[i]*np.conjugate(cube_b.cube[i]))/np.sqrt(energies_a[i]*energies_b[i]))
    corr_err = np.sqrt(1 - corr**2)/ np.sqrt(nx*ny/cube_a.scales**2 - 2)
    
    obj_name,filt_a,filt_b = give_names(cube_a.name,cube_b.name)
    plt.clf()
    plt.errorbar(cube_a.scales,corr,yerr = corr_err,fmt = '.')
    plt.xlabel(r'Scale (pixels)', fontsize=14)
    plt.ylabel(r'Correlation', fontsize=14)
    plt.title(r'Scale Correlation of '+obj_name+' in '+filt_a+' and '+filt_b+' Filters', fontsize=16)
    plt.savefig(path+'\\'+'Output/'+cube_a.name+'_'+cube_b.name+'_pethat_corr_smooth.png', dpi=400)
    plt.show()
    return