#	REMOVE COSMIC RAY WITH ASTROSCRAPPY
#	https://github.com/astropy/astroscrappy/blob/master/README.rst
#	https://astroscrappy.readthedocs.io/en/latest/
#	conda install -c conda-forge astroscrappy 
#------------------------------------------------------------
#	CREATED BY CHANGSU CHOI
#	2020.01.21	MODIFIED BY Gregory S.H. Paek
#============================================================
import astropy.io.fits as fits
import astroscrappy as cr
import os, glob
from imsng import tool
#============================================================
#	FUNCTION
#============================================================
def crclean(inim, outim, gain=1.0, rdnoise=6.5, sigclip=4.5, sigfrac=0.3, objlim=5.0, niter=4, cleantype='medmask', fsmode='median', verbose=True):
	data, hdr = fits.getdata(inim, header=True)
	c1,c2=cr.detect_cosmics(data, gain=gain, readnoise=rdnoise,
							sigclip=sigclip, sigfrac=sigfrac, objlim=objlim, niter=niter,
							cleantype=cleantype, fsmode=fsmode,
							# psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765,
							verbose=verbose)

	hdr.set('COMMENT','AstroSCRAPPY cr removed')
	fits.writeto(outim, c2, header=hdr)



path_obs = '/home/sonic/Research/table'

# os.system("ls Calibrated*.fits > obj.list")
# os.system('rm cCal*.fits')
os.system('ls *.fits *.fit')
imlist = glob.glob(input('IMAGES TO PROCESS\t: ')); imlist.sort()


for inim in imlist:
	try:
		obs = inim.split('-')[1]
	except:
		obs = input('CCD?\t: ')
	gain, pixscale, rdnoise = tool.getccdinfo(obs, path_obs)
	param_crclean = dict(inim = inim,
						outim = 'c'+inim,
						gain = gain, rdnoise=rdnoise,
						sigfrac=0.3, objlim=5.0, niter=4, cleantype='medmask', fsmode='median',
						verbose=True)
	crclean(**param_crclean)


# for im in lists : 
	# ascrapp(im)
	# os.system('ls cCalibrated*.fits | wc -l && ls Calibrated*.fits |wc -l') 


# print 'done','\a'