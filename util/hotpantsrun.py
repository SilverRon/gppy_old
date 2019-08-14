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
import alipy

#============================================================
def hotpants(imlist, refim):
    if type(imlist) != list: imlist = [imlist]
    starttime	= time.time()
    for inim in imlist:
        outfile = os.path.dirname(inim)+'/hd'+os.path.basename(inim)
        convfile= os.path.dirname(inim)+'/hc'+os.path.basename(inim)
        #com     = 'hotpants -c t -n i -iu 60000 -tu 60000 -v 0 -inim '+infile[n]+' -tmplim '+ref_img+' -outim '+outfile+' -oci '+convfile
        #com='hotpants -v 0 -c i -n i -inim '+infile[n]+' -tmplim ref.fits -outim '+outfile+' -oci '+convfile
        #com='hotpants -c t -n i -inim '+infile[n]+' -tmplim ref.fits -outim '+outfile+' -oci '+convfile
        #com     = 'hotpants -c t -n i -iu 60000 -tl -40yu0 -tu 1000000 -v 0 -inim '+infile[n]+' -tmplim '+ref_img+' -outim '+outfile+' -oci '+convfile
        com     = 'hotpants -c t -n i -iu 60000 -tu 60000 -tl -10000 -v 0 -inim '+inim+' -tmplim '+refim+' -outim '+outfile+' -oci '+convfile
        os.system(com)
    deltime		= time.time() - starttime
    print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
#------------------------------------------------------------
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
