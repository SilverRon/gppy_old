from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.table import Table

photbl = ascii.read('phot.dat')
#------------------------------------------------------------
bandlist = []
for inim in photbl['obs']:
	bandlist.append(inim.split('-')[5])
photbl['band'] = np.array(bandlist)

#------------------------------------------------------------
objlist = []
for inim in photbl['obs']:
	objlist.append(inim.split('-')[2])
photbl['object'] = np.array(objlist)

datelist	= []
for inim in photbl['obs']:
	part	= inim.split('-')
	datelist.append('20'+part[3][0:2]+'-'+part[3][2:4]+'-'+part[3][4:6]+'T'+part[4][0:2]+':'+part[4][2:4]+':'+part[4][4:])
photbl['date-obs'] = np.array(datelist)

jdlist		= []
for date in photbl['date-obs']:
	jdlist.append(Time(date, format='isot', scale='utc').jd)
photbl['jd'] = np.array(jdlist)

#photbl['obs'] = 'UKIRT'

newtbl		= Table()
newtbl['object'], newtbl['obs'], newtbl['band'], newtbl['date-obs'], newtbl['jd'], newtbl['seeing'], newtbl['mag'] = photbl['object'], photbl['obs'], photbl['band'], photbl['date-obs'], photbl['jd'], photbl['seeing'], photbl['mag']
newtbl.write('phot_ukirt.dat', format='ascii', overwrite=True)
#------------------------------------------------------------





photbl[ (photbl['band']=='R') & (photbl['seeing']<3.0) & (photbl['mag']>20.0)]
photbl[ (photbl['band']=='V') & (photbl['seeing']<4.0) & (photbl['mag']>19.0) ]
photbl[ (photbl['band']=='B') & (photbl['seeing']<2.5) & (photbl['mag']>20.6) ]

import os, glob
path_qso	= '/data3/IMSNG/LOAO/gal'
path_gundam	= '/mnt/window/Users/User/Downloads/data/loao/ref'
downcom		= 'sshpass -prjseka23! scp -ro StrictHostKeyChecking=no paek@qso.snu.ac.kr:'

#reflist		= ['NGC0488', 'NGC1309', 'UGC02855', 'NGC2207', 'NGC2993', 'IC2537', 'NGC3169', 'NGC3294', 'NGC3344', 'NGC3629', 'NGC3646', 'NGC3938', 'NGC4030', 'NGC4108']
#reflist			= ['NGC3169', 'NGC3294', 'NGC3344', 'NGC3629', 'NGC3646', 'NGC3938']
#reflist			= ['NGC0488', 'NGC1309', 'NGC2207', 'NGC2993', 'UGC02855', 'IC2537', 'NGC1385']
reflist		= ['M95', 'M66', 'M98', 'M86', 'M87', 'M91', 'M90', 'M58', 'M64', 'M63']

for obj in reflist:
	os.system('mkdir '+path_gundam+'/'+obj+'/')
	com		= downcom+path_qso+'/'+obj+'/C*.fits '+path_gundam+'/'+obj+'/'
	print(com)
	os.system(com)