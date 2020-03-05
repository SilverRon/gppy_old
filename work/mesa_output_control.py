#	MESA OUTPUT CONTROL
#	CREATED	2019.11.04	Gregory S.H. Paek
#============================================================
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from imsng import phot
import time
# import mesa_reader to make its classes accessible
import mesa_reader as mr
#============================================================
#	FUNCTION
#============================================================

#============================================================
#	USER SETTING
#============================================================
path_base	= '/home/sonic/CourseWork/19B/StellarEV/term_project/test'
#------------------------------------------------------------
# make a MesaData object from a history file
h = mr.MesaData(path_base+'/LOGS/history.data')

# extract the star_age column of data
ages = h.data('star_age')
# nages = ages/np.max(ages)/1e6		#	[Myr]

plt.close('all')
plt.scatter(h.log_Teff, h.log_L, c=ages)
plt.xlabel(r'$\log T_{eff}$', fontsize=20)
plt.ylabel(r'$\log L$', fontsize=20)
plt.colorbar()

plt.minorticks_on()
plt.tick_params(which='both', direction='in')
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.gca().invert_xaxis()
plt.tight_layout()

plt.savefig(path_base+'/plot/hr.png', overwrite=True)

'''
l = mr.MesaLogDir(path_base+'/LOGS')


# load the profile associated with model number 100
p_100 = l.profile_data(100)
# the same as the following
p_100 = l.profile_data(model_number=100)

# load the profile with PROFILE number 12
p_12 = l.profile_data(profile_number=12)

# load the last profile saved (largest model number)
p_last = l.profile_data()
'''
#------------------------------------------------------------
#	IMAGES TO PHOTOMETRY
#------------------------------------------------------------

#============================================================
#	MAIN COMMAND
#============================================================

#------------------------------------------------------------
#	FINISH
#------------------------------------------------------------