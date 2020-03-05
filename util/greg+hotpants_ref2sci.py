#============================================================
#	PYTHON SCRIPT FOR SUBTRACTION USING HOTPANTS
#	Usage : 
#	OR JUST RUN AND WRITE INPUT & REF IMAGE
#	2019.07.05	GREGORY S.H. PAEK
#============================================================
import numpy as np
import os, sys, glob
from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u
from multiprocessing import Process, Pool
import multiprocessing as mp
import time
import alipy
#============================================================
def hotpants(imlist, refim, ngmode=False, insig=0, refsig=0):
	if type(imlist) != list: imlist = [imlist]
	starttime	= time.time()
	for inim in imlist:
		outfile = os.path.dirname(inim)+'/hd'+os.path.basename(inim)
		convfile= os.path.dirname(inim)+'/hc'+os.path.basename(inim)
		if ngmode == False:
			com     = 'hotpants -c t -n i -iu 60000 -tu 60000 -tl -10000 -v 0 -inim '+inim+' -tmplim '+refim+' -outim '+outfile+' -oci '+convfile
		elif ngmode == True:
			sigmatch	= np.sqrt(insig**2-refsig**2)
			com		= 'hotpants -c t -n i -iu 60000 -tu 60000 -tl -10000 -v 0 -inim {0} -tmplim {1} -outim {2} -oci {3} -ng 3 6 {4} 4 {5} 2 {6}'.format(inim, refim, outfile, convfile, 0.5*sigmatch, sigmatch, 2.0*sigmatch)
		os.system(com)
	deltime		= time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
#-------------------------------------------------------------
def gregistering(images_to_align, ref_image):
	starttime	= time.time()
	if type(images_to_align) != list: images_to_align = [images_to_align]
	if ref_image == '': ref_image = images_to_align[0]
	identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
	for id in identifications: # list of the same length as images_to_align.
		if id.ok == True: # i.e., if it worked
			print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
		else:
			print "%20s : no transformation found !" % (id.ukn.name)
	outputshape = alipy.align.shape(ref_image)
	for id in identifications:
		if id.ok == True:
			params_align	= dict(	filepath	= id.ukn.filepath,
									uknstarlist	= id.uknmatchstars,
									refstarlist	= id.refmatchstars,
									shape		= alipy.align.shape(ref_image),
									outdir		= './',
									makepng		= False)
			alipy.align.irafalign(**params_align)
	deltime		= time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
#-------------------------------------------------------------
def subt_routine(inim, refim, ngmode=False, insig=0, refsig=0):
	#time.sleep(5)
	grefim	= inim[:-5]+'_ref2sci.fits'
	gregistering(refim, inim)
	os.system('mv {0} {1}'.format(refim[:-5]+'_gregister.fits', grefim))
	hotpants(inim, grefim, ngmode, insig, refsig)
#============================================================
path_base	= '.'
refim	= path_base+'/ref.fits'
imlist	= glob.glob(path_base+'/Calib-*com.fits')	;imlist.sort()
# imlist	= glob.glob(path_base+'/Calib-*0.fits')	;imlist.sort()
#------------------------------------------------------------
for inim in imlist:
	# insig	= fits.getheader(inim)['seeing']/2.355
	# refsig	= fits.getheader(refim)['seeing']/2.355
	subt_routine(inim, refim)#, ngmode=True, insig=insig, refsig=refsig)
#------------------------------------------------------------
#	MULTI PROCESSING
#------------------------------------------------------------
'''
for inim in imlist:
	if __name__ == '__main__':
		jobs	= []
		p		= mp.Process(target=subt_routine, args=(inim, refim))
		jobs.append(p)
		p.start()
		p.join()

if __name__ == '__main__':
	procs	= []
	for inim in imlist:
		proc = mp.Process(target=subt_routine, args=(inim, refim, ))
		procs.append(proc)
		proc.start()
	for proc in procs:
		proc.join()
'''