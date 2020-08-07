#	ABSOLUTE MAG -> APPARENT MAG WITH DISTANCE INFO.
#============================================================
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from astropy.table import Table, vstack
from astropy import units as u

def abs2app(M, Mer, d, der):
	m = M + 5 * np.log10(d) - 5
	mer = np.sqrt( (Mer)**2 + ((5*der)/(d*np.log(10))) )
	return m, mer



path_table = '/home/sonic/Research/gppy/table/phot_gw170817_Drout.abs.dat'
path_save = '/home/sonic/Research/gppy/table'

photbl = ascii.read(path_table)

dist, dister = 156 * u.Mpc, 41 * u.Mpc
d, der = dist.to(u.pc).value, dister.to(u.pc).value
M, Mer = photbl['mag_abs'], photbl['mager_abs']

m, mer = abs2app(M, Mer, d, der)




photbl['mag_gw190425'] = m
photbl['mager_gw190425'] = mer

photbl.write('{}/phot_gw170817_Drout.gw190425.dat'.format(path_save), format='ascii.tab', overwrite=True)