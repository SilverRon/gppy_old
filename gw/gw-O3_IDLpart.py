#   GW ALERT FOR O3 RUN IN QSO SERVER
#   BASED ON https://emfollow.docs.ligo.org/userguide/tutorial/index.html
#	PLEASE RUN DIRECTORY /data3/gwshare/bayestar/IDLzone
#   2019.03.03	MADE	BY Gregory S.H. Paek
#	2019.03.28	UPDATED	BY Gregory S.H. Paek
#	2019.04.18	UPDATED	BY Gregory S.H. Paek
#	2019.04.27	UPDATED	BY Gregory S.H. Paek
#	2019.09.25 REPAIRED	BY Gregory S.H. Paek
#	2019.10.30	UPDATED	BY Gregory S.H. Paek
#	2019.11.11	UPDATED BY Gregory S.H. Paek
#============================================================
#	MODULE
#------------------------------------------------------------
import gcn
import gcn.handlers
import gcn.notice_types
import healpy as hp
import numpy as np
import warnings
warnings.filterwarnings(action='ignore')
print('{0}\nIDL PART FOR GW ALERT IS STANDBY.\n{1}\n\n\n'.format('='*60, '='*60))
'''
import lxml.etree
path_test = '/home/gw/Research/test'
payload		= open(path_test+'/MS181101ab-1-Preliminary.xml', 'rb').read()
root		= lxml.etree.fromstring(payload)
'''
#============================================================
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
	#------------------------------------------------------------
	#	Read all of the VOEvent parameters from the "What" section.
	#------------------------------------------------------------
	starttime	= time.time()
	params	= {elem.attrib['name']:
			  elem.attrib['value']
			  for elem in root.iterfind('.//Param')}
	role	= root.attrib['role']
	#------------------------------------------------------------
	eventname	= params['GraceID']+'-'+params['AlertType']
	path_idlzone	= '/data3/gwshare/bayestar/IDLzone'
	# path_healpix	= path_idlzone+'/'+eventname+'.fits.gz'
	path_healpix	= path_idlzone+'/bayestar.fits.gz'
	path_save		= '/data3/gwshare/bayestar/'+eventname
	os.system('mkdir '+path_save)
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
		progenitor	= 'UNKNOWN'
		pass
	#------------------------------------------------------------
	#	IDL PART (Joonho Kim)
	#------------------------------------------------------------
	'''
	outlist = [	'gwfollowup_kmtnet_target_{}.txt'.format(eventname.split('-')[0]),
				'gwfollowup_kmtnet_target_{}_30.txt'.format(eventname.split('-')[0]),
				'sc{}_B.cat'.format(eventname.split('-')[0]),
				'sc{}_R.cat'.format(eventname.split('-')[0]),
				'sc{}_30R.cat'.format(eventname.split('-')[0])]
	'''
	# if root.attrib['role'] != 'test':
	if root.attrib['role'] != 'test':
		print('# DOWNLOAD BAYESTER FILE')
		healpixfits	= params['skymap_fits']
		os.system('wget {0} -O {1}'.format(healpixfits, path_healpix))
		#------------------------------------------------------------
		os.system('gzip -d {}'.format(path_healpix))
		os.system('chmod 777 {}'.format(path_healpix[:-3]))
		f = open(path_idlzone+'/gwinfo.dat', 'w')
		f.write('{}\n{}\n'.format(eventname, path_save))
		f.close()
		#------------------------------------------------------------
		#	IDL BODY
		#------------------------------------------------------------
		os.system('sh /data3/gwshare/idlpro/gwfollowup_kmtnet.sh > /data3/gwshare/idlpro/gwfollowup_kmtnet.sh.log 2>&1')
		print('CLEANING & MOVING DUMMY FILES')
		# os.system('mv {}/{} {}'.format(path_idlzone, eventname+'-kmtnet*.txt', path_save))
		os.system('rm {}/gwinfo.dat {}/a.out {}/bayestar.fits'.format(path_idlzone, path_idlzone, path_idlzone))
		os.system('chmod 777 {}/*'.format(path_save))
		'''
		os.system('rm {}'.format(path_idlzone+'/*.fits'))
		os.system('rm {}'.format(path_idlzone+'/gwinfo.dat'))
		for output in outlist:
			print('mv {} {}'.format(path_idlzone+'/'+output, path_save+'/'))
			os.system('mv {} {}'.format(path_idlzone+'/'+output, path_save+'/'))
		'''
	else:
		print('TEST SIGNAL')
	deltime		= time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)'+'\n')
	print('IDL PART FOR GW ALERT IS STANDBY.')
#------------------------------------------------------------
#   Listen for GCNs until the program is interrupted
#   (killed or interrupted with control-C).
gcn.listen(handler=process_gcn)
#------------------------------------------------------------