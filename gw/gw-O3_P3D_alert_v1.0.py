#   UTILIZE GW P_3D LOCALIZATION MAP
#   BASED ON https://github.com/lpsinger/gw-galaxies/blob/master/gw-galaxies.ipynb
#   2019.07.19	MADE BY		Gregory S.H. Paek
#	2019.XX.XX	UPDATED	BY	Gregory S.H. Paek
#	2019.08.28	MODIFIED BY Gregory S.H. Paek
#	2019.10.30	UPDATED BY	Gregory S.H. Paek
#	2019.11.11	UPDATED BY	Gregory S.H. Paek
#============================================================
#	MODULE
#------------------------------------------------------------
import gcn
import gcn.handlers
import gcn.notice_types
import healpy as hp
import numpy as np
import time
import os, glob
from astropy.table import Table, Column, MaskedColumn, vstack
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import matplotlib.pyplot as plt
from imsng import gw
from imsng import tool
from ligo.skymap.postprocess import find_greedy_credible_levels
#------------------------------------------------------------
import astropy.utils.data
#------------------------------------------------------------
#	FOR THE GALAXY CROSS MATCHING
from astropy.table import Table, vstack, hstack, Column
import ligo.skymap.plot
from scipy.stats import norm
import scipy.stats
import warnings
warnings.filterwarnings(action='ignore')
#============================================================
print('{0}\nMAIN PART FOR GW ALERT IS STANDBY.\n{1}\n\n\n'.format('='*60, '='*60))

#	TEST ALERT
'''
import lxml.etree
path_test = '/home/gw/Research/test'
payload		= open(path_test+'/MS181101ab-1-Preliminary.xml', 'rb').read()
root		= lxml.etree.fromstring(payload)
'''
@gcn.handlers.include_notice_types(
	gcn.notice_types.LVC_PRELIMINARY,
	gcn.notice_types.LVC_INITIAL,
	gcn.notice_types.LVC_UPDATE)
def process_gcn(payload, root):
	import gcn
	import gcn.handlers
	import gcn.notice_types
	import healpy as hp
	import numpy as np
	import time
	import os, glob, sys
	from astropy.table import Table, Column, MaskedColumn, vstack
	import astropy.units as u
	from astropy.coordinates import SkyCoord
	from astropy.io import ascii
	import matplotlib.pyplot as plt
	from imsng import gw
	from imsng import tool
	from ligo.skymap.postprocess import find_greedy_credible_levels
	#------------------------------------------------------------
	import astropy.utils.data
	#------------------------------------------------------------
	#	FOR THE GALAXY CROSS MATCHING
	from astropy.table import Table, vstack, hstack, Column
	import ligo.skymap.plot
	from scipy.stats import norm
	import scipy.stats
	import warnings
	warnings.filterwarnings(action='ignore')
	'''
	EXTRACT CONFIDENCE MAP FROM BAYESTAR FILE
	'''
	#------------------------------------------------------------
	#	INITIAL SETTING
	#------------------------------------------------------------
	starttime = time.time()
	yymmdd, hhmmss = tool.timename()
	yy = yymmdd[0:2]
	yyyy = '20'+yy
	mm, dd = yymmdd[2:4], yymmdd[4:6]
	#------------------------------------------------------------
	#	PATH
	#------------------------------------------------------------
	path_catalog = '/home/gw/Research/GLADE_2.3+2MASS_PSC+identi_name.dat'
	path_gw = '/data3/gwshare/bayestar'
	#------------------------------------------------------------
	#	INITIAL PARAMETERS
	#------------------------------------------------------------
	# confidence = 0.9
	# conflist = [0.5, 0.7, 0.9, 0.99]
	alldeg = (4*np.pi)*((180/np.pi)**2)
	distsigcut = 3
	xmlkeys = ['GraceID', 'AlertType', 'EventPage', 'Instruments', 'FAR', 'Group', 'skymap_fits', 'BNS', 'NSBH', 'BBH', 'Terrestrial', 'HasNS', 'HasRemnant']
	hdrkeys = [ 'NSIDE',
				'OBJECT',
				'REFERENC',
				'INSTRUME',
				'DATE-OBS',
				'MJD-OBS',
				'DATE',
				'ORIGIN',
				'RUNTIME',
				'DISTMEAN',
				'DISTSTD',
				'LOGBCI',
				'LOGBSN',
				'DATE-BLD',
				'HISTORY']
	#------------------------------------------------------------

	# print('{0}\n{1}\n{2}'.format('='*60, '0-0.\tGW ALERT', '-'*60))
	#------------------------------------------------------------
	#	Read all of the VOEvent parameters from the "What" section.
	#------------------------------------------------------------
	params = {	elem.attrib['name']:
				elem.attrib['value']
				for elem in root.iterfind('.//Param')}
	role = root.attrib['role']
	eventname = params['GraceID']+'-'+params['AlertType']
	path_healpix = path_gw+'/'+eventname+'/'+eventname+'.fits.gz'
	#------------------------------------------------------------
	progs = np.array(['BNS', 'NSBH', 'BBH', 'MassGap', 'Terrestrial'])
	progs_prob = []
	for key in progs:
		if key in params.keys():
			progs_prob.append(float(params[key]))
			indx_prog = np.where(progs_prob == np.max(progs_prob))
		else:
			pass
	prog = np.asscalar(progs[indx_prog])
	#------------------------------------------------------------
	path_save = path_gw + '/' + eventname
	os.system('mkdir '+path_save)
	sys.stdout = open(path_save+'/log.txt', 'a')
	print('{0}\n{1}\tGW ALERT\n{2}'.format('='*60, yymmdd+'-'+hhmmss, '-'*60))
	print('# EVENTNAME :\t{0} ({1}) by {2}'.format(eventname, role, prog))
	#------------------------------------------------------------
	#	DOWNLOADE BAYESTER FITS
	#------------------------------------------------------------
	print('# DOWNLOAD BAYESTER FILE')
	healpixfits	= params['skymap_fits']
	os.system('wget {0} -O {1}'.format(healpixfits, path_healpix))
	#	PLOT HEALpix SKYMAP
	plotcom = 'ligo-skymap-plot {} -o {} --annotate --contour 50 90'.format(path_healpix, path_save+'/'+eventname+'.png')
	os.system(plotcom)
	try:
		#------------------------------------------------------------
		#	1. bayestar.fits -> RA/Dec/P_2D TABLE
		#------------------------------------------------------------
		cantbl, prob, confareainfo, hdr = gw.read3Dhealpix2candidates(path_healpix, path_catalog, eventname=eventname, distsigcut=distsigcut, header=True)
		npix = len(prob)
		cantbl['sort'] = 'K'
		cantbl['rank'] = np.arange(len(cantbl))+1
		cantbl['cumsum'] = np.cumsum(cantbl['Prob'])
		cantbl.write(path_save+'/'+eventname+'-all_candi.dat', format='ascii', overwrite=True)
		#------------------------------------------------------------#
		#	PLOT
		#------------------------------------------------------------#
		gw.plotcumscore(cantbl, probcutlist=[0.5, 0.90, 0.95, 0.99], eventname=eventname, path_save=path_save)
		gw.expectedLC(eventname, hdr, path_save)
		#------------------------------------------------------------#
		#	MAKE REGION FILE
		#------------------------------------------------------------#
		param_regmaker	= dict(	filename = path_save+'/'+eventname+'.reg',
								name = cantbl['ra'],
								ra = cantbl['dec'],
								dec = cantbl['iden_name'])
		tool.ds9regmaker(**param_regmaker)
		print('# MAKE REGION FILE (.reg)')
		#------------------------------------------------------------#
		#	WRITE SUMMARY
		#------------------------------------------------------------#
		print('# MAKE SUMMARY FILE (.txt)')
		f	= open(path_save+'/'+eventname+'-summary.txt', 'w')
		f.write('Role :\t'+root.attrib['role']+'\n')
		for key in xmlkeys:
			f.write(key+' :\t'+str(params[key])+'\n')
			if key == 'FAR':
				f.write('1/FAR :\t1/{} yr\n'.format(int(1/float(params['FAR'])/(365*24*60*60))))
		f.write('='*60+'\n')
		for key in hdrkeys:
			f.write(key+' :\t'+str(hdr[key])+'\n')
		for conf in [0.5, 0.9]:
			param_heal	= dict(	healpixfits=healpixfits,
								confidence=conf,
								eventname=eventname,
								save_path=path_save,
								hdr=True,
								view=False)
			healtbl, hdr= gw.heal2table(**param_heal)
			skydeg_cut	= (alldeg/npix)*len(healtbl)
			skydeg_cut_N= (alldeg/npix)*len(healtbl[healtbl['dec']>0])
			skydeg_cut_S= (alldeg/npix)*len(healtbl[healtbl['dec']<0])
			f.write(str(conf*100)+'% conf. area :\t '+str(round(skydeg_cut, 1))+'\t[deg2]\n')
			f.write(str(conf*100)+'% conf. area (N.sphere) :\t '+str(round(skydeg_cut_N, 1))+'\t[deg2]\n')
			f.write(str(conf*100)+'% conf. area (S.sphere) :\t '+str(round(skydeg_cut_S, 1))+'\t[deg2]\n')
		indx_max = np.where(healtbl['P_2D'] == np.max(healtbl['P_2D']))
		indx_min = np.where(healtbl['P_2D'] == np.min(healtbl['P_2D']))
		cmax = SkyCoord(healtbl['ra'][indx_max], healtbl['dec'][indx_max], frame='icrs', unit='deg')
		cmaxstr = (cmax.to_string('hmsdms')[0]).split(' ')
		ramax, demax = cmaxstr[0], cmaxstr[1]
		# cmin = SkyCoord(healtbl['ra'][indx_min], healtbl['dec'][indx_min], frame='icrs', unit='deg')
		# cminstr = (cmin.to_string('hmsdms')[0]).split(' ')
		# ramin, demin = cminstr[0], cminstr[1]
		f.write('Number of candidates :\t{}\n'.format(len(cantbl)))
		for conf in [0.99, 0.9, 0.5]:
			f.write('Number of candidates (cumscore<{}):\t{}\n'.format(conf, len(cantbl[cantbl['cumsum']<=conf])))
		mag, magerr	= gw.kilonova_mag(hdr['DISTMEAN'], hdr['DISTSTD'])
		f.write('EXPECTED KN APP.MAG IN i-BAND [AB] : '+str(round(mag, 2))+'+/-'+str(round(magerr, 2))+'\n')
		f.write('MAX_RA,DEC :\t{} {}\n'.format(ramax, demax))
		# f.write('MIN_RA,DEC :\t{} {}\n'.format(ramin, demin))
		#	TILE NUMBER FOR KMTNet
		try:
			tilefile = open(path_save+'/'+eventname+'-kmtnet.txt')
			tilenumb = len(tilefile.readlines())-4
			f.write('KMTNet TILE # :\t{}\n'.format(tilenumb))
		except:
			print('FAIL TO READ IDL RESULTS')
		# f.close()
		os.system('chmod 777 '+path_save+'/*')
		#------------------------------------------------------------#
		#	6.	MAKE TARGET LIST FOR EACH OBSERVATORIES
		#------------------------------------------------------------#
		print('# MAKE TARGET LIST FOR EACH OBSERVATORIES (.txt)')
		if params['AlertType'] == 'Preliminary':
			deldays	= 0
			tblcut = 0.50
		elif params['AlertType'] == 'Initial':
			deldays = 1
			tblcut = 0.95
		elif params['AlertType'] == 'Update':
			deldays = 2
			tblcut = 0.99
		elif role == 'test':
			deldays = 0
			tblcut = 0.50
		else:
			deldays = 0
			tblcut = 0.50

		cutbl = cantbl[cantbl['cumsum']<tblcut]
		rtstbl = Table()
		rtstbl['name'], rtstbl['ra'], rtstbl['dec'], rtstbl['sort'], rtstbl['rank'], rtstbl['dist']	= cutbl['iden_name'], cutbl['ra'], cutbl['dec'], cutbl['sort'], cutbl['rank'], cutbl['dist']
		rtstbl.meta['cumscorecut'] = tblcut
		rtstbl.write(path_save+'/'+eventname+'-rts_candi_in{}.txt'.format(tblcut), format='ascii', overwrite=True)
		obspath = '/home/gw/Research/observatory.txt'
		path_catalog = path_save+'/'+eventname+'-rts_candi_in{}.txt'.format(tblcut)

		# for obs in ['SAO', 'McD', 'MAO', 'LOAO', 'SSO', 'UKIRT', 'Spain', 'NMex', 'Calif']:
		# for obs in ['SAO', 'McD', 'MAO', 'LOAO', 'SSO', 'UKIRT']:
		for obs in ['SAO', 'McD', 'LOAO', 'SSO', 'UKIRT']:
			# print(obs)
			param_rts	= dict(	observatory=obs,
								headname=eventname,
								save_path=path_save,
								obspath=obspath,
								catpath=path_catalog,
								start=yyyy+'/'+mm+'/'+dd,
								end=yyyy+'/'+mm+'/'+str(int(dd)+deldays),
								altlimit=30.,
								moonseperation=40.,
								sunlimit='-18',
								numlimit=300)
			tool.rtsmaker(**param_rts)
			rtsfiles = glob.glob('{}/{}-*-rts_vis-{}.txt'.format(path_save, eventname, obs))
			for rts in rtsfiles:
				rtsfile = open(rts)
				if len(rtsfile.readlines()) <= 11:
					obsout = '{}\tIS NOT OBSERVABLE LOCATION.'.format(obs)
					print(obsout)
					os.system('rm {}'.format(rts))
				else:
					obsout = '{}\tIS OBSERVABLE LOCATION.'.format(obs)
					print(obsout)
					pass
				f.write(obsout+'\n')
		f.close()
	except:
		print('THERE IS A ERROR DURING PROCESS.\nSEND MAIL TO STAFFS.')
	#------------------------------------------------------------#
	#	6.5 EXCEPTION FOR SQUEAN (REMOVE ':')
	#------------------------------------------------------------#
	'''
	rtsfile	= path_save+'/'+eventname+'-'+yyyy+mm+dd+'-rts_vis-McD.txt'
	obsinfo	= ascii.read(rtsfile, header_start=10, data_start=11)
	obsname	= (obsinfo['name'])
	obsra	= (obsinfo['ra'])
	obsdec	= (obsinfo['dec'])

	f		= open(directory+"radec_"+stryear+strmonth+strday+"_"+observatory.lower()+".txt","w")
	for n in range(len(obsname)):
		obsra[n]	= obsra[n].replace(':',' ')
		obsdec[n]	= obsdec[n].replace(':',' ')
		f.write('{:s} {:s} {:s}'.format(obsname[n],obsra[n],obsdec[n])+'\n')
	f.close()
	'''
	#------------------------------------------------------------#
	#	7.	SEND MAIL
	#------------------------------------------------------------#
	if (role != 'test') & (float(params['Terrestrial']) < 0.7):
		# print('{0}\n# SEND MAIL TO CO-WORKERS {1}'.format('-'*60, len(cantbl), '-'*60))
		print('{0}\n# SEND MAIL TO CO-WORKERS\n{1}'.format('-'*60, '-'*60))
		reciver	= ascii.read('/home/gw/Research/alert-mail-reciver.txt')
		subject	= '[GW alert-{}] {} @{}+/-{}Mpc'.format(prog, eventname, int(hdr['DISTMEAN']), int(hdr['DISTSTD']))
		# contents= 'CONTENTS'
		import codecs
		try:
			contents= codecs.open(path_save+'/'+eventname+'-summary.txt', 'rb', 'utf-8')
		except:
			contents='NONE'
		fromID	= 'ceouobs@gmail.com'
		fromPW	= 'ceou@snu'
		toIDs	= ''
		for address in reciver['address']: toIDs += address+','
		# toIDs	= toIDs[:-1]
		# toIDs	= "gregorypaek94@gmail.com"
		# toIDs	= 'joonhoyo305@gmail.com'
		ccIDs	= 'gundam_psh@naver.com'
		import glob
		# path	= glob.glob(path_save+'/'+eventname+'-*.txt')
		path	= []
		for typ in ('*.txt', '*.cat'):
			typfiles = glob.glob('{}/*{}-{}'.format(path_save, eventname, typ))
			path.extend(typfiles)
		#path.append(path_healpix_tmp+yymmdd+'/'+yymmdd+'-'+hhmmss+'-skymap.png')
		#tool.send_gmail(subject, fromID, contents, fromPW, toIDs, ccIDs=None, path=path)
		import os
		import smtplib
		from email.mime.base import MIMEBase
		from email.mime.text import MIMEText
		from email.mime.image import MIMEImage
		from email.mime.multipart import MIMEMultipart
		from email.header import Header  
		#msg		= MIMEBase('mixed')
		#msg		= MIMEText(contents, 'plain', 'utf-8')
		msg		= MIMEMultipart()
		msg['Subject']	= Header(s=subject, charset="utf-8")
		msg['From']		= fromID
		msg['To']		= toIDs
		if ccIDs != None:
			msg['Cc']		= ccIDs
		#msg.attach(MIMEText(contents, 'plain', 'utf-8'))
		msg.attach(MIMEText(contents.read()))

		#	ATTACH TEXT FILE ON MAIL
		if path != None:
			if type(path) != list:
				filelist	= []
				filelist.append(path)
			else:
				filelist	= path
			for file in filelist:
				part	= MIMEBase("application", "octet-stream")
				part.set_payload(open(file, 'rb').read())
				part.add_header(	'Content-Disposition',
									'attachment; filename="%s"'% os.path.basename(file))
				msg.attach(part)
		#	ATTACH PNG FILE ON MAIL
		pnglist		= glob.glob(path_save+'/'+eventname+'*.png')
		for png in pnglist:
			fp		= open(png, 'rb')
			img		= MIMEImage(fp.read())
			fp.close()
			img.add_header('Content-Disposition', 'attachment', filename=os.path.basename(png))
			msg.attach(img)
		#	ACCESS TO GMAIL & SEND MAIL
		smtp_gmail		= smtplib.SMTP_SSL('smtp.gmail.com', 465)
		smtp_gmail.login(fromID, fromPW)
		smtp_gmail.sendmail(msg["From"], msg["To"].split(",") + msg["Cc"].split(","), msg.as_string())
		smtp_gmail.quit()
		comment	= 'Send '+str(path)+'\nFrom\t'+fromID+'\nTo'; print(comment); print(toIDs)
	# os.system('cat '+path_save+'/'+eventname+'-summary.txt')
	os.system('chmod 777 '+path_save)
	os.system('chmod 777 '+path_save+'/*')
	deltime		= time.time() - starttime
	print('DONE ({0} sec)'.format(round(deltime, 1)))
	print('MAIN PART FOR GW ALERT IS STANDBY.\n{0}\n\n\n'.format('='*60))
#------------------------------------------------------------#
#   Listen for GCNs until the program is interrupted
#   (killed or interrupted with control-C).
gcn.listen(handler=process_gcn)
#------------------------------------------------------------#
