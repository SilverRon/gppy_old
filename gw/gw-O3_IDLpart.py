#   GW ALERT FOR O3 RUN IN QSO SERVER
#   BASED ON https://emfollow.docs.ligo.org/userguide/tutorial/index.html
#	PLEASE RUN DIRECTORY /data3/gwshare/bayestar/IDLzone
#   2019.03.03 MADE		BY Gregory S.H. Paek
#	2019.03.28 UPDATED	BY Gregory S.H. Paek
#	2019.04.18 UPDATED	BY Gregory S.H. Paek
#	2019.04.27 UPDATED	BY Gregory S.H. Paek
#============================================================#
#	MODULE
#------------------------------------------------------------#
import gcn
import gcn.handlers
import gcn.notice_types
import healpy as hp
import numpy as np
#============================================================#
# Function to call every time a GCN is received.
# Run only for notices of type
# LVC_PRELIMINARY, LVC_INITIAL, or LVC_UPDATE.
@gcn.handlers.include_notice_types(
	gcn.notice_types.LVC_PRELIMINARY,
	gcn.notice_types.LVC_INITIAL,
	gcn.notice_types.LVC_UPDATE)
def process_gcn(payload, root):
	'''
	EXTRACT CONFIDENCE MAP FROM BAYESTAR FILE
	'''
	import gcn
	import gcn.handlers
	import gcn.notice_types
	import healpy as hp
	import numpy as np
	import time
	import os, glob
	from astropy.table import Table, Column, MaskedColumn
	import astropy.units as u
	from astropy.coordinates import SkyCoord
	from astropy.io import ascii
	import matplotlib.pyplot as plt
	from imsng import gw
	from imsng import tool
	from astropy.table import Table, vstack
	#------------------------------------------------------------#
	#	Read all of the VOEvent parameters from the "What" section.
	#------------------------------------------------------------#
	#	GIVE DELAY TO ESCAPE OVERAPPED COMMAND
	#time.sleep(10)
	starttime	= time.time()
	params	= {elem.attrib['name']:
			  elem.attrib['value']
			  for elem in root.iterfind('.//Param')}
	role	= root.attrib['role']
	#	Print all parameters.
	eventname	= params['GraceID']+'_'+params['AlertType']
	#save_healfix	= '/data3/gwshare/bayestar/'
	#save_path		= save_healfix+eventname
	save_healfix	= '/data3/gwshare/bayestar/IDLzone'
	save_path		= '/data3/gwshare/bayestar/'+eventname

	os.system('mkdir '+save_path)

	yymmdd, hhmmss	= tool.timename()
	yy	= yymmdd[0:2]
	yyyy= '20'+yy
	mm	= yymmdd[2:4]
	dd	= yymmdd[4:6]
	#------------------------------------------------------------#
	#	PROGENITOR
	try:
		prog		= ['BNS', 'NSBH', 'BBH', 'Terrestrial']
		prog_prob	= []
		for key in prog:
			value	= params[key]
			prog_prob.append(float(value))
		indx_prog	= np.where(np.max(prog_prob) == prog_prob)
		progenitor	= prog[indx_prog[0][0]]
	except:
		pass
	#------------------------------------------------------------#
	healpixfits	= params['skymap_fits']
	os.system('wget '+healpixfits)

	#	IDL PART (Joonho Kim)
	if root.attrib['role'] != 'test':
		gzlist	= glob.glob(save_healfix+'/*.gz')
		os.system('gzip -d '+gzlist[0])
		os.system('sh /data3/gwshare/idlpro/gwfollowup_kmtnet.sh')
		baylist	= glob.glob(save_healfix+'/*.fits')
		os.system('mv '+baylist[0]+' '+save_path+'/'+eventname+'-bayestar.fits')
		os.system('mv '+save_healfix+'/sc*.cat '+save_path+'/')
		os.system('mv '+save_healfix+'/gw*.txt '+save_path+'/')
		os.system('chmod 777 '+save_path+'/*')
	else:
		print('TEST SIGNAL')
		os.system('rm bayestar.fits.gz')
	deltime		= time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
	print('IDL PART FOR GW ALERT IS STANDBY.')

'''
#	TEST ALERT
import lxml.etree
payload		= open('/home/gw/Research/test/MS181101ab-1-Preliminary.xml', 'rb').read()
root		= lxml.etree.fromstring(payload)

#	TEST BAYESTAR
confidence	= 0.9
healpixfits	= '/home/sonic/Research/cat/GWlocal/bayestar170814-refined.fits'
progenitor	= 'BBH'
'''

#------------------------------------------------------------#
#   Listen for GCNs until the program is interrupted
#   (killed or interrupted with control-C).
gcn.listen(handler=process_gcn)
#------------------------------------------------------------#
