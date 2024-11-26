import numpy as np
import scipy, time
import matplotlib
import pathlib
import os

from alive_progress import alive_bar
import matplotlib.animation as animation

from scipy import integrate, fft ,signal
import pylab as plt
from scipy import io
from astropy.io import fits
from scipy import stats
from typing import Literal

#-----------------------------------------------------------------
#---------------------- Pre-Processing ---------------------------
#-----------------------------------------------------------------

path = str(pathlib.Path(__file__).parent.resolve())

def save_array_as_FITS(array,filename):
    hdu_new = fits.PrimaryHDU(array)
    hdu_new.writeto(path+'\\data\\'+filename+'.fits',overwrite=True)
    return

def preprocess(name,filename,crop_cor = []):
    # Retrieves observational data
    outname = 'Output/'+name+'/'+name+'_'
    adata = fits.open(path+'\\data\\'+filename)
    a = adata[0].data
    
    # Crop the image if it's required
    if crop_cor != []:
        a=a[crop_cor[0][0]:crop_cor[0][1],crop_cor[1][0]:crop_cor[1][1]]

    kont = np.array(a.astype(float))
    kont = kont/np.max(kont)  # Normalize to 1
    og = kont

    # Creates a new folder if necessary
    newpath = path+'\\'+ 'Output/'+ name
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    
    # Plot and save the image of the object
    plt.clf()
    plt.imshow(kont, origin='lower', interpolation='nearest')
    plt.title("Normalized "+name)
    plt.colorbar()
    plt.savefig(path+"\\"+outname +'original.png', dpi=200)
    #plt.show()
    plt.clf()

    # Zero-padding the image to a square to reduce the complications with Pet Hat later on
    ny = kont[0,:].size
    nx = kont[:,0].size
    kont = np.pad(kont,[(0,int(np.max([-nx+ny,0]))), (0,int(np.max([nx-ny,0])))], mode='median') 

    return kont, nx, ny, og # Returns the final map and the shape of the original image

#-----------------------------------------------------------------
#---------------------- The Pet-Hat function ---------------------
#-----------------------------------------------------------------

def pethat_k(k): # This function is in the fourier space 
    if ((np.pi <= k)and(k <= (4.0*np.pi))): 
        return np.cos(np.pi * 0.5 * np.log2(k * 0.5/np.pi))**2
    else:
        return 0.


def pethat_phi_func(data,a_scale): # Builds an array of a 2d Pet Hat map using scale
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


def pethat_wavelet_scale_analysis(name,filename,crop_cor=[],scales_in=[2,100,1], scales_type = "triplet", pixel_scale = 1, distance = 1 , full = False):  # placeholder for pixel scale and distance
    scales_types = ["triplet", "array"]
    if scales_type not in scales_types:
        raise ValueError("Invalid scale type. Expected one of: %s" % scales_types)
    begin_str = "Commencing the wavelet scale analysis of "+name
    end_str = "Wavelet scale analysis of "+name+" has finished."

    print(begin_str)
    kont, nx, ny, og = preprocess(name,filename,crop_cor)

    if full:
        # Produces a complete FFT of the image
        fft_shape = 2*np.array(kont.shape) - 1
        fft_kont = scipy.fft.fftshift(scipy.fft.fft2(kont,fft_shape))

    if full == False:
        # Produces a complete FFT of the image
        fft_kont = scipy.fft.fftshift(scipy.fft.fft2(kont))
    
    # Calculates the coefficient needed to convert the unit of the sacle to pixels
    fnx = fft_kont[:,0].size
    coeff = 4*np.pi/fnx 

    # Defines the array for scales and wavelet coefficients
    if scales_type == "triplet":
        [start,n,step_length] = scales_in
        en_scales = start + np.arange(0,n) * step_length
    if scales_type == "array":
        en_scales = scales_in
        n = len(en_scales)
    wavelet_coeffs = np.clongdouble(np.zeros((n,nx,ny)))

    # Calculates the wavelet coefficients for all the scales
    with alive_bar(n) as bar: 
        for i in range(n):
            a = en_scales[i] * coeff # The scale
            map1 = scipy.fft.ifft2(scipy.fft.ifftshift(fft_kont*pethat_phi_func(fft_kont,a)))
            map1 = map1[0:nx,0:ny]
            wavelet_coeffs[i] = map1
            bar()
    print(end_str)
    return wavelet_results(wavelet_coeffs,en_scales,name,og,pixel_scale,distance)

#-----------------------------------------------------------------
#------------------ Wavelet Results as an Object -----------------
#-----------------------------------------------------------------

      
class wavelet_results:
    """ 
    This object consists of the name of the observed object, the cube of wavelet coefficients, all the corresponding scales, 
    the original image, the pixel scale in arcsec/pixel, and the distance to the object in parsecs. The object can then
    produce all the needed data, utilizing several methods. 

    The possible (at least for now) actions consist of:
    creating a gif from the data cube of wavelet coeffs, 
    saving the data cube in a FITS file, 
    saving the wavelet coeffs as PNG images, 
    calculating the energies and plot the energy over scale plot, 
    and calculating the sum of the wavelet coeffs. 
    """

    def __init__(self, cube, scales , name , og, pixel_scale, dist):
        self.cube = cube
        self.scales = scales
        self.name = name
        self.outname = 'Output/'+name+'/'+name+'_'
        self.original = og  
        self.pixel_scale = pixel_scale
        self.dist = dist
    
    def unit_conv(self, to_unit = 'pixels'):
        if to_unit == 'pixels':
            return 1
        if to_unit == 'arcsec':
            return self.pixel_scale
        if to_unit == 'arcmin':
            return self.pixel_scale/60
        if to_unit == 'pc':
            return self.dist*self.pixel_scale/206265
        if to_unit == 'kpc':
            return self.dist*self.pixel_scale/206265/1000
        else:
            raise ValueError("Invalid unit. Expected one of: %s" % ['pixels','arcsec','arcmin','pc','kpc'])
        
    def create_gif(self,unit='pixels'):
        def update(frame,cube,scales,name):
            img.set_array(cube[frame])
            tx.set_text('PetHat wavelet animation of ' + name + ', scale = '+str(np.round(scales[frame]*self.unit_conv(unit),3))+' '+ unit)
        
        print("Creating GIF from the coefficients of "+self.name+" ... ",end="")

        our_data = np.real(self.cube)
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(1,1,1)
        img = ax.imshow(our_data[0], origin='lower', interpolation='nearest',vmin=np.min(our_data), vmax=np.max(our_data), animated=True)
        cb = fig.colorbar(img)
        tx = ax.set_title('PetHat wavelet animation of ' + self.name + ', scale = '+str(np.round(self.scales[0]*self.unit_conv(unit),3))+' '+ unit, fontsize=12)

        ani = animation.FuncAnimation(fig=fig, func=update, fargs=(our_data,self.scales,self.name), frames=len(self.cube)-1, interval=100)
        ani.save(filename=path+"\\"+self.outname +'animation.gif', writer='ffmpeg',codec="libx264")
        print("DONE")
        return
    
    def save_FITS(self):
        print("Saving the data of "+self.name+" as a FITS file ... ",end="")
        hdu_new = fits.PrimaryHDU()
        hdu_ima = fits.ImageHDU(np.real(self.cube))
        hdu_ens = fits.ImageHDU(self.calc_energies(do_plot=False))

        # Update header
        hdu_ima.header.append(('AUTHOR','Aria Hemmatian'))
        hdu_ima.header.append(('OBJECT',self.name))
        hdu_ima.header.append(('SCA_UNIT','pixels', 'unit of the scales of the layers'))
        hdu_ima.header.append(('DIST',self.dist,'distance of the object in parsecs'))
        hdu_ima.header.append(('PIX_SCA',self.pixel_scale,'pixel scale in arcsecs'))

        hdul = fits.HDUList([hdu_new,hdu_ima,hdu_ens])
        hdul.writeto(path+'\\'+self.outname+'wavelet_coeffs.fits',overwrite=True)
        print("DONE")
        return
        
    def save_layers(self,unit='pixels'):
        print("Saving the layers of "+self.name+" as PNG ... ",end="")
        for i in range(len(self.cube)):
            plt.imshow(np.real(self.cube[i]), origin='lower', interpolation='nearest')
            plt.colorbar()
            plt.title('PetHat wavelet coefficient of '+self.name+', scale = '+str(np.round(self.scales[i]*self.unit_conv(unit),3))+' '+ unit, fontsize=12)
            plt.savefig(path+"\\"+self.outname +'wavelet_pethat_fast_'+str(i+1)+'.png', dpi=300)
            plt.clf()
        print("DONE")
        return

    def calc_energies(self, do_plot=True, unit = 'pixels'):
        energies = np.zeros(len(self.cube))
        for i in range(len(self.cube)):
            energies[i] = np.sum(np.abs(self.cube[i])**2)
        if do_plot:
            plt.clf()
            plt.plot(self.scales*self.unit_conv(unit),energies)
            plt.xlabel(r'Scale ('+unit+')', fontsize=11)
            plt.ylabel(r'Wavelet Energy of '+self.name, fontsize=11)
            plt.title('PetHat Wavelet Energies', fontsize=12)
            plt.savefig(path+'\\'+self.outname+'pethat_energy_smooth.png', dpi=300)
            #plt.show()
            plt.clf()
        return energies
    
    def calc_sum(self, do_plot=True):
        sums = np.sum(np.real(self.cube),axis=0)
        if do_plot:
            plt.imshow(sums, origin='lower', interpolation='nearest')
            plt.colorbar()
            plt.title('Sum of the PetHat Wavelet Coefficients', fontsize=12)
            plt.savefig(path+"\\"+self.outname +'wavelet_pethat_fast_total.png', dpi=200)
            #plt.show()
            plt.clf()
        return sums


#-----------------------------------------------------------------
#----------------- Wavelet Scale Cross-Correlation ---------------
#-----------------------------------------------------------------

def give_names(str_a,str_b): # This function extracts the name of the object and the filters for the cross-correlation function
    name = ''
    for i in range(len(str_a)):
        if str_a[i] == '_':
            str_a = str_a[i+1:]
            str_b = str_b[i+1:]
            break
        name+=(str_a[i])  
    return name , str_a, str_b

def plot_correlation(cube_a, cube_b, unit='pixels'): # This function gets two result objects and if they are compatible, plots their scale cross-correlation over scale
    energies_a = cube_a.calc_energies(do_plot=False)
    energies_b = cube_b.calc_energies(do_plot=False)
    if (not np.allclose(cube_a.scales,cube_b.scales)) or np.shape(cube_a.cube) != np.shape(cube_b.cube):
        print('Error: Cubes are not compatible.')
        return
    obj_name,filt_a,filt_b = give_names(cube_a.name,cube_b.name)
    print("The normal correlation between the maps is: ", (np.sum( cube_a.original * np.conjugate(cube_b.original) )/np.sqrt( np.sum(np.abs(cube_a.original)**2) * np.sum(np.abs(cube_b.original)**2) )))
    (n,nx,ny) = np.shape(cube_a.cube)
    
    # Calculating the cross-correlation and its error
    corr = np.zeros(n)
    for i in range(n):
        corr[i] = np.real(np.sum(cube_a.cube[i]*np.conjugate(cube_b.cube[i]))/np.sqrt(energies_a[i]*energies_b[i]))
    corr_err = np.sqrt(1 - corr**2)/ np.sqrt(nx*ny/cube_a.scales**2 - 2)
    
    plt.clf()
    plt.errorbar(cube_a.scales*cube_a.unit_conv(unit),corr,yerr = corr_err,fmt = '.')
    plt.xlabel(r'Scale ('+unit+')', fontsize=12)
    plt.ylabel(r'Correlation', fontsize=12)
    plt.title(r'Scale Correlation of '+obj_name+' in '+filt_a+' and '+filt_b+' Filters', fontsize=12)
    plt.savefig(path+'\\'+'Output/'+obj_name+'_'+filt_a+'_'+filt_b+'_pethat_corr_smooth.png', dpi=400)
    #plt.show()
    plt.clf()
    return