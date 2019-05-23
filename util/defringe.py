#	FRINGE PROCESS BASED ON PROF. IM's CODE
#	WORK FOR ONLY LOAO DATA (OR YOU NEED TO MAKE 'fringe_i_ori.fits' and fringe_i.dat)
#	19.03.12	GREGORY S.H. PAEK
#	19.05.21	GREGORY S.H. PAEK
#============================================================
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import astropy.io.fits as fits
#============================================================
def fringe_cal(inim):
	data, hdr	= fits.getdata(inim, header=True)
	#subframe	= data[512:1537, 512:1537]
	#skymean, skymed, skysig		= bkgest_mask(inim)
	#subdata		= data-skymed

	dfr_list	= []
	for n in range(len(fringes)):
		'''
		fringe_b= np.median(subdata[yb1[n]:yb2[n], xb1[n]:xb2[n]])
		fringe_f= np.median(subdata[yf1[n]:yf2[n], xf1[n]:xf2[n]])
		'''
		fringe_b= np.median(data[yb1[n]:yb2[n], xb1[n]:xb2[n]])
		fringe_f= np.median(data[yf1[n]:yf2[n], xf1[n]:xf2[n]])		
		dfringe	= fringe_b-fringe_f
		dfr_list.append(dfringe)
	return np.array(dfr_list)


path_base	= '/home/sonic/Research/yourpy/fringe'

masterf	= path_base+'/fringe_i_ori.fits'
dataf, hdrf	= fits.getdata(masterf, header=True)

imlist	= glob.glob('Calib*.fits')
imlist.sort()
fringes	= ascii.read(path_base+'/fringe_i.dat')

xb1, xb2= fringes['xb']-5, fringes['xb']+5
yb1, yb2= fringes['yb']-5, fringes['yb']+5

xf1, xf2= fringes['xf']-5, fringes['xf']+5
yf1, yf2= fringes['yf']-5, fringes['yf']+5


master_fri		= fringe_cal(masterf)
infr_list		= []
for inim in imlist:
	fscale		= np.median(fringe_cal(inim)/master_fri)
	#print(fscale)
	data, hdr	= fits.getdata(inim, header=True)
	fri_scaled	= dataf*fscale
	newim		= data-fri_scaled
	fits.writeto('Df'+inim, newim, hdr, overwrite=True)
	#print(np.median(data), np.median(fri_scaled), np.median(newim))