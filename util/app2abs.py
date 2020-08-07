#	APPARENT MAG -> ABSOLUTE MAG WITH DISTANCE INFO.
#============================================================
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from astropy.table import Table, vstack
from astropy import units as u

def app2abs(m, mer, d, der):
	M = m -5 * np.log10(d) + 5
	Mer = np.sqrt( (mer)**2 + ((-5*der)/(d*np.log(10)))**2 )
	return M, Mer



path_table = '/home/sonic/Research/gppy/table/phot_gw170817_Drout.dat'
path_save = '/home/sonic/Research/gppy/table'





photbl = ascii.read(path_table)

dist, dister = 37.7 * u.Mpc, 8.7 * u.Mpc
d, der = dist.to(u.pc).value, dister.to(u.pc).value
m, mer = photbl['mag'], photbl['magerr']

M = m -5 * np.log10(d) + 5
Mer = np.sqrt( (mer)**2 + ((-5*der)/(d*np.log(10)))**2 )

photbl['mag_abs'] = M
photbl['mager_abs'] = Mer

photbl.write('{}/phot_gw170817_Drout.abs.dat'.format(path_save), format='ascii.tab', overwrite=True)