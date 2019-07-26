#	PHOTOMETRY CODE FOR PYTHON 3.X, AND BULK OF KMTNet IMAGES
#	CREATED	2019.06.20	Gregory S.H. Paek
#	MODIFIED 2019.07.26 Gregory S.H. Paek
#============================================================
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.time import Time
#from multiprocessing import Process, Pool
#import multiprocessing as mp
import time
#============================================================
photbl = ascii.read('phot.dat')
iphtbl, uphtbl = photbl, photbl
initbl = ascii.read('inthefov_initial.dat')
updtbl = ascii.read('inthefov_update.dat')
#------------------------------------------------------------
sourini, sourupd = [], []
for inim in photbl['image']:
	indxini = np.where(iphtbl['image']==inim)
	sourini.append(np.copy(initbl['sources'][indxini])[0])
	indxupd = np.where(uphtbl['image']==inim)
	sourupd.append(np.copy(updtbl['sources'][indxupd])[0])
#------------------------------------------------------------
sourini, sourupd = np.array(sourini), np.array(sourupd)
iphtbl['sources'], uphtbl['sources'] = sourini, sourupd
#------------------------------------------------------------
iphtbl.write('phot_ini_sources.dat', format='ascii')
uphtbl.write('phot_upd_sources.dat', format='ascii')
