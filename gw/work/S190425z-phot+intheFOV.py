#	PHOTOMETRY CODE FOR PYTHON 3.X, AND BULK OF KMTNet IMAGES
#	CREATED	2019.06.20	Gregory S.H. Paek
#	MODIFIED 2019.07.26 Gregory S.H. Paek
#============================================================
import os
import glob
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.time import Time
import time
photbl = ascii.read('phot.dat')
#============================================================
#	INITIAL
#------------------------------------------------------------
iphtbl = photbl
initbl = ascii.read('inthefov_initial.dat')
#------------------------------------------------------------
sourini = []
for inim in photbl['image']:
	indxini = np.where(iphtbl['image']==inim)
	sourini.append(np.asscalar(initbl['sources'][indxini]))
iphtbl['sources'] = np.array(sourini)
#------------------------------------------------------------
if 'phot_ini_sources.dat' in glob.glob('*.dat'):
	os.system('mv phot_ini_sources.dat phot_ini_sources.dat.bkg')
iphtbl.write('phot_ini_sources.dat', format='ascii')
#============================================================
#	UPDATE
#------------------------------------------------------------
uphtbl = photbl
updtbl = ascii.read('inthefov_update.dat')
#------------------------------------------------------------
sourupd = []
for inim in photbl['image']:
	indxupd = np.where(uphtbl['image']==inim)
	sourupd.append(np.asscalar(updtbl['sources'][indxupd]))
uphtbl['sources'] = np.array(sourupd)
#------------------------------------------------------------
if 'phot_upd_sources.dat' in glob.glob('*.dat'):
	os.system('mv phot_upd_sources.dat phot_upd_sources.dat.bkg')
uphtbl.write('phot_upd_sources.dat', format='ascii')









#------------------------------------------------------------
#	COMBINE ALL TABLES
#------------------------------------------------------------
#	INITIAL
tblist = []
for tbl in glob.glob('phot*sources*.dat'):
    tblist.append(ascii.read(tbl))
comtbl = vstack(tblist)
comtbl.write('phot_ini_sources_GECKO.dat', format='ascii')
#	UPDATE
tblist = []
for tbl in glob.glob('phot*sources*.dat'):
	tblist.append(ascii.read(tbl))
comtbl = vstack(tblist)
comtbl.write('phot_upd_sources_GECKO.dat', format='ascii')
