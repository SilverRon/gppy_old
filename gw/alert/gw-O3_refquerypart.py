#   GW ALERT FOR ER 14 RUN IN QSO SERVER
#   BASED ON https://emfollow.docs.ligo.org/userguide/tutorial/index.html
#   2019.03.03 MADE		BY Gregory S.H. Paek
#	2019.03.28 UPDATED	BY Gregory S.H. Paek
#	2019.04.18 UPDATED	BY Gregory S.H. Paek
#	2019.05.05 UPDATED	BY Gregory S.H. Paek
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
	from astropy.io import ascii
	import matplotlib.pyplot as plt
	from imsng import gw
	from imsng import tool
	from imsng import query
	#------------------------------------------------------------#
	time.sleep(180)
	starttime	= time.time()
	confidence	= 0.9
	#------------------------------------------------------------#
	#	Read all of the VOEvent parameters from the "What" section.
	#------------------------------------------------------------#
	params	= {elem.attrib['name']:
			  elem.attrib['value']
			  for elem in root.iterfind('.//Param')}

	role	= root.attrib['role']

	if role != 'test':
		eventname	= params['GraceID']+'_'+params['AlertType']

		save_healfix	= '/data3/gwshare/bayestar/'
		save_path		= save_healfix+eventname
		save_PS1		= '/data3/gwshare/bayestar/PS1'
		targettbl		= ascii.read(save_path+'/'+eventname+'-all_candi.txt')

		existlist		= glob.glob(save_PS1+'/Ref-PS1*.fits')

		fail_list		= []
		#for i in range(len(targettbl)):
		for i in range(300):
			print('['+str(i+1)+'/'+str(len(targettbl))+']')
			if 'Ref-PS1-'+targettbl['name'][i]+'-r.fits' not in existlist:
				param_down= dict(	name		= targettbl['name'][i],
									ra			= targettbl['ra'][i],
									dec			= targettbl['dec'][i],
									size		= 5000,
									output_size	= None,
									filters		= 'r',
									format		= "fits",
									save_dir	= save_PS1)
				try:
					query.downimage_routine(**param_down)
				except:
					try:
						query.downimage_routine(**param_down)
					except:
						try:
							query.downimage_routine(**param_down)
						except:
							print('NOT COVERED BY PS1.')
							fail_list.append(targettbl['name'][i])
			else:
				pass
		f	= open(save_path+'/'+eventname+'-query_from_PS1_fail.txt', 'w')
		f.write('name\n')
		for fail in fail_list:
			f.write(fail+'\n')	
		f.close()

	deltime		= time.time() - starttime
	print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
	print('PS1 QUERY PART FOR GW ALERT IS STANDBY.')

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
