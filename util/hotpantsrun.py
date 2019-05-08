#============================================================
#	PYTHON SCRIPT FOR SUBTRACTION USING HOTPANTS
#	Usage : 
#	python hotpantsrun.py str1 ref_image
#	OR JUST RUN AND WRITE INPUT & REF IMAGE
#	2019.05.08	GREGORY S.H. PAEK
#============================================================
import numpy as np
import os, sys, glob
from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u
from multiprocessing import Process, Pool
import multiprocessing as mp
import time
#============================================================
def hotpants(imlist, refim):
    starttime	= time.time()
    for inim in imlist:
        outfile = 'hd' + inim
        convfile= 'hc' + inim
        #com     = 'hotpants -c t -n i -iu 60000 -tu 60000 -v 0 -inim '+infile[n]+' -tmplim '+ref_img+' -outim '+outfile+' -oci '+convfile
        #com='hotpants -v 0 -c i -n i -inim '+infile[n]+' -tmplim ref.fits -outim '+outfile+' -oci '+convfile
        #com='hotpants -c t -n i -inim '+infile[n]+' -tmplim ref.fits -outim '+outfile+' -oci '+convfile
        #com     = 'hotpants -c t -n i -iu 60000 -tl -40yu0 -tu 1000000 -v 0 -inim '+infile[n]+' -tmplim '+ref_img+' -outim '+outfile+' -oci '+convfile
        com     = 'hotpants -c t -n i -iu 60000 -tu 60000 -tl -10000 -v 0 -inim '+inim+' -tmplim '+refim+' -outim '+outfile+' -oci '+convfile
        os.system(com)
    deltime		= time.time() - starttime
    print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
#------------------------------------------------------------
#	INPUT
#------------------------------------------------------------
os.system('ls *.fits')
try     :
    if sys.argv[1] == 'on'  : command0 = 'Cal*ter.fits'
except  : command0    = raw_input('Image list to process (Cal*ter.fits): ')
if command0 == '' : command0    = 'Cal*ter.fits'

try     :
    if sys.argv[1] == 'on'  : command1    = 'ref.fits'
except  : command1    = raw_input('Reference Image (ref.fits) : ')
if command1 == '' : command1    = 'ref.fits'
refim       = command1
imlist      = glob.glob(command0)
#------------------------------------------------------------
#	MULTI PROCESSING
#------------------------------------------------------------
if __name__ == '__main__':
	jobs	= []
	p			= mp.Process(target=hotpants, args=(imlist, refim))
	jobs.append(p)
	p.start()