#   GW ALERT FOR ER 14 RUN IN QSO SERVER
#   BASED ON https://emfollow.docs.ligo.org/userguide/tutorial/index.html
#   2019.03.03 MADE		BY Gregory S.H. Paek
#	2019.03.28 UPDATED	BY Gregory S.H. Paek
#	2019.04.18 UPDATED	BY Gregory S.H. Paek
#	2019.05.02 UPDATED	BY Gregory S.H. Paek (ver. 2.0)
#	2019.05.10 UPDATED	BY Gregory S.H. Paek
#	2019.07.08 UPDATED	BY Gregory S.H. Paek (ver. 3.0)
#============================================================#
#	MODULE
#------------------------------------------------------------#
import gcn
import gcn.handlers
import gcn.notice_types
import healpy as hp
import numpy as np
#============================================================#
print('{0}\nMAIN PART FOR GW ALERT IS STANDBY.\n{1}\n\n\n'.format('='*60, '='*60))
'''
#	TEST ALERT
import lxml.etree
payload		= open('/home/gw/Research/test/MS181101ab-1-Preliminary.xml', 'rb').read()
root		= lxml.etree.fromstring(payload)
'''
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
	import os
	from astropy.table import Table, Column, MaskedColumn
	import astropy.units as u
	from astropy.coordinates import SkyCoord
	from astropy.io import ascii
	import matplotlib.pyplot as plt
	from imsng import gw
	from imsng import tool
	from astropy.table import Table, vstack
	#------------------------------------------------------------
	#	INITIAL SETTING
	#------------------------------------------------------------
	starttime		= time.time()
	yymmdd, hhmmss	= tool.timename()
	yy				= yymmdd[0:2]
	yyyy			= '20'+yy
	mm, dd			= yymmdd[2:4], yymmdd[4:6]
	#------------------------------------------------------------
	catpath			= '/home/gw/Research/GLADE_2.3+2MASS_PSC+identi_name.dat'
	save_healfix	= '/data3/gwshare/bayestar/'
	confidence		= 0.5
	print('{0}\n{1}\n{2}'.format('='*60, '0-0.\tGW ALERT', '-'*60))
	#------------------------------------------------------------#
	#	Read all of the VOEvent parameters from the "What" section.
	#------------------------------------------------------------#
	params	= {	elem.attrib['name']:
				elem.attrib['value']
				for elem in root.iterfind('.//Param')}
	role	= root.attrib['role']
	eventname	= params['GraceID']+'-'+params['AlertType']
	#------------------------------------------------------------#
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
	print('0-1.\tEVENTNAME :\t{0} ({1}) by {2}'.format(eventname, role, progenitor))
	save_path		= save_healfix+eventname
	os.system('mkdir '+save_path)
	#------------------------------------------------------------#
	#	DOWNLOADE BAYESTER FITS
	#------------------------------------------------------------#
	print('0-2.\tDOWNLOAD BAYESTER FILE')
	healpixfits	= params['skymap_fits']
	os.system('wget {0}'.format(healpixfits))
	os.system('mv {0} {1}'.format(os.path.basename(healpixfits), save_path+'/'+eventname+'.fits.gz'))
	#------------------------------------------------------------#
	#	1. bayestar.fits -> RA/Dec/P_2D TABLE
	#------------------------------------------------------------#
	print('{0}\n1-0.\tHOST GALAXY CANDIDATES\n{1}'.format('-'*60, '-'*60))
	print('1-1.\tREAD BAYESTER FILE')
	param_heal	= dict(	healpixfits=healpixfits,
						confidence=0.9,
						eventname=eventname,
						save_path=save_path,
						hdr=True,
						view=True)
	healtbl, hdr= gw.heal2table(**param_heal)
	healtbl.write(save_path+'/'+eventname+'_healpix.dat', format='ascii', overwrite=True)
	npix        = hdr['NAXIS2']
	gwdist      = hdr['DISTMEAN']
	gwdiststd   = hdr['DISTSTD']
	updist, lodist	= gwdist+2*gwdiststd, gwdist-2*gwdiststd
	if lodist < 0: lodist = 0
	#------------------------------------------------------------#
	#	READ GLADE catalog
	#------------------------------------------------------------#
	print('1-2.\tREAD GLADE CATALOG')
	cattbl		= ascii.read(catpath)
	cattbl		= cattbl[ (cattbl['dist']<updist) & (cattbl['dist']>lodist) ]
	#------------------------------------------------------------#
	#	3.	MATCHING GLADE & HEALPIX
	#------------------------------------------------------------#
	print('1-3.\tMATCH WITH GLADE')
	#	SKY COORD.
	coo_intbl   = SkyCoord(healtbl['ra'], healtbl['dec'], unit=(u.deg, u.deg))
	coo_reftbl  = SkyCoord(cattbl['ra'], cattbl['dec'], unit=(u.deg, u.deg))
	#   INDEX FOR REF.TABLE
	indx_match, d2d, d3d  = coo_intbl.match_to_catalog_sky(coo_reftbl, nthneighbor=1)
	#	REMOVE MATCHING & null & dist.cut
	alldeg		= (4*np.pi)*((180./np.pi)**2)		#	[deg2]
	maxerr	    = 2.*np.sqrt(alldeg / npix)*3600.	#	["/pix]
	dummytbl		= cattbl[indx_match]
	dummytbl['sep']	= d2d
	dummytbl['P_2D']= healtbl['P_2D']
	matchtbl0		= dummytbl[	(dummytbl['sep']*3600. < maxerr) ]
	#	REMOVE OVERAP
	voidlst			= []
	for i in matchtbl0['dist']:
		voidlst.append(i)
	distlst			= list(set(voidlst))
	indx_min		= []
	for i in distlst:
		indx_dum	= np.where(i == matchtbl0['dist'])
		indx_alive	= np.where(np.min(matchtbl0['sep'][indx_dum]) == matchtbl0['sep'])
		indx_min.append(indx_alive[0][0])
	matchtbl		= matchtbl0[indx_min]
	#------------------------------------------------------------#
	#	SEPERATE TABLE
	#------------------------------------------------------------#
	print('{0}\n2-0.\tSCORE {1} CANDIDATES \n{2}'.format('-'*60, len(matchtbl), '-'*60))
	#	1. K TABLE / 2. NO K BUT B TABLE / 3. NO K AND B TABLE (ONLY DIST)
	cantblK			= matchtbl[	(matchtbl['K']!=-99)]
	cantblB			= matchtbl[	(matchtbl['K']==-99) &
								(matchtbl['B']!=-99)]
	cantbldist		= matchtbl[	(matchtbl['K']==-99) &
								(matchtbl['B']==-99)]
	print('2-1.\tAVAILABLE : K {0} / NO K BUT B {1} / ONLY DIST {2}'.format(len(cantblK), len(cantblB), len(cantbldist)))
	#------------------------------------------------------------#
	#	SORT: P_2D & Luminosity & Distance
	#------------------------------------------------------------#
	#	PARAMETERS
	param_K			= dict(	gwdist=gwdist, 
							gwdister=gwdiststd,
							gldist=cantblK['dist'],
							gldister=cantblK['dist_err'],
							appmag=cantblK['K'],
							prob=cantblK['P_2D'])
	param_B			= dict(	gwdist=gwdist, 
							gwdister=gwdiststd,
							gldist=cantblB['dist'],
							gldister=cantblB['dist_err'],
							appmag=cantblB['B'],
							prob=cantblB['P_2D'])
	param_dist		= dict(	gwdist=gwdist, 
							gwdister=gwdiststd,
							gldist=cantbldist['dist'],
							gldister=cantbldist['dist_err'],
							appmag=np.ones(len(cantbldist)),
							prob=cantbldist['P_2D'])
	#------------------------------------------------------------#
	try:
		tbllist	= []
		if	len(cantblK) > 0:
			# print('K-band info. is available.')
			cantblK['score']	= gw.prob(**param_K)
			cantblK				= cantblK[np.argsort(-1*cantblK['score'])]
			cantblK['sort']		= 'K'
			tbllist.append(cantblK)
		if	len(cantblB) > 0:
			# print('B-band info. is available.')
			cantblB['score']	= gw.prob(**param_B)
			cantblB				= cantblB[np.argsort(-1*cantblB['score'])]
			cantblB['score']	= cantblB['score']*( np.min(cantblK['score'])*0.99/np.max(cantblB['score']) )
			cantblB['sort']		= 'B'
			tbllist.append(cantblB)
		if	len(cantbldist) > 0:
			# print('Distance info. is available.')
			cantbldist['score']	= gw.prob(**param_dist)
			cantbldist			= cantbldist[np.argsort(-1*cantbldist['score'])]
			try:
				cantbldist['score']	= cantbldist['score']*( np.min(cantblB['score'])*0.99/np.max(cantbldist['score']) )
				cantbldist['sort']	= 'dist'
				tbllist.append(cantbldist)
			except:
				cantbldist['score']	= cantbldist['score']*( np.min(cantblK['score'])*0.99/np.max(cantbldist['score']) )
				cantbldist['sort']	= 'dist'
				tbllist.append(cantbldist)
		cantbl				= vstack(tbllist)
		cantbl['score']		= cantbl['score']/np.sum(cantbl['score'])
	except:
		print('\n**FAIL TO SORTING PROCESS. MAKE TABLE WITHOUT CONSTRAINTS.**\n')
		cantbl		= Table(	[matchtbl['name'], matchtbl['ra'], matchtbl['dec']], matchtbl['P_2D']/np.sum(matchtbl['P_2D']),
								names=['name', 'ra', 'dec', 'score'])
		cantbl['sort']= 'none'
	cantbl['rank']		= np.arange(len(cantbl))+1
	cantbl.write(save_path+'/'+eventname+'-all_candi.txt', format='ascii', overwrite=True)
	targtbl			= Table(	[cantbl['name'], cantbl['ra'], cantbl['dec'], cantbl['sort']],
									names=['name', 'ra', 'dec', 'sort'])
	print('2-2.\tALL SCORED CANDIDATES {0}'.format(len(targtbl)))
	#------------------------------------------------------------#
	#	PLOT
	#------------------------------------------------------------#
	gw.expectedLC(eventname, hdr, save_path)
	plt.close('all')
	cumscore	= np.cumsum(cantbl['score'])
	plt.plot(	np.arange(len(cantbl)), cumscore, color='tomato',
				label=str(len(cantbl))+' targets')
	plt.axhline(y=0.9, color='dodgerblue', linestyle=':',
				label=str(0.9)+' ('+str(len(cumscore[cumscore<0.9]))+')')
	plt.axhline(y=0.5, color='green', linestyle=':',
				label=str(0.5)+' ('+str(len(cumscore[cumscore<0.5]))+')')
	plt.title(eventname+' '+str(confidence*100)+'%')
	plt.legend()
	plt.savefig(save_path+'/'+eventname+'-cumscore.png', dpi=500)
	print('2-3.\tPLOT CUMULATIVE SCORE & EXPECTED LC (.png)')
	#------------------------------------------------------------#
	#	MAKE REGION FILE
	#------------------------------------------------------------#
	param_regmaker	= dict(	filename	= save_path+'/'+eventname,
							intbl		= targtbl,
							name		= 'name',
							racol		= 'ra',
							deccol		= 'dec')
	tool.ds9regmaker(**param_regmaker)
	print('2-4.\tMAKE REGION FILE (.reg)')
	#------------------------------------------------------------#
	#	WRITE SUMMARY
	#------------------------------------------------------------#
	print('2-5.\tMAKE SUMMARY FILE (.txt)')
	f	= open(save_path+'/'+eventname+'-summary.txt', 'w')
	f.write('\nRole\t: '+root.attrib['role']+'\n')
	for key, value in params.items():
		f.write(key+'\t'+str(value)+'\n')
	f.write('#'+'-'*60+'#'+'\n')
	for key in hdr.keys():
		line = key+'\t'+str(hdr[key])+'\n'
		f.write(line)
	for conf in [0.5, 0.9]:
		param_heal	= dict(	healpixfits=healpixfits,
							confidence=conf,
							eventname=eventname,
							save_path=save_path,
							hdr=True,
							view=True)
		healtbl, hdr= gw.heal2table(**param_heal)
		skydeg_cut	= (alldeg/npix)*len(healtbl)
		skydeg_cut_N= (alldeg/npix)*len(healtbl[healtbl['dec']>0])
		skydeg_cut_S= (alldeg/npix)*len(healtbl[healtbl['dec']<0])
		f.write(str(conf*100)+'% confidence area\t: '+str(skydeg_cut)+'\t[deg2]\n')
		f.write(str(conf*100)+'% confidence area for Northern sphere\t: '+str(skydeg_cut_N)+'\t[deg2]\n')
		f.write(str(conf*100)+'% confidence area for Southern sphere\t: '+str(skydeg_cut_S)+'\t[deg2]\n')

	mag, magerr	= gw.kilonova_mag(gwdist, gwdiststd)
	f.write('EXPECTED KILONOVA APP.MAG IN i-BAND [AB] : '+str(round(mag, 2))+'+/-'+str(round(magerr, 2))+'\n')
	f.write('Number of candidates (90% confidence)\t:{0}\n'.format(len(matchtbl)))
	f.write('CENTER_RA\t: '+str(healtbl[ healtbl['P_2D'] == np.max(healtbl['P_2D']) ]['ra'][0])+'\n')
	f.write('CENTER_Dec\t: '+str(healtbl[ healtbl['P_2D'] == np.max(healtbl['P_2D']) ]['dec'][0])+'\n')
	f.write('MAX_RA\t: '+str(np.max(healtbl['ra']))+'\n')
	f.write('MAX_Dec\t: '+str(np.max(healtbl['dec']))+'\n')
	f.write('MIN_RA\t: '+str(np.min(healtbl['ra']))+'\n')
	f.write('MIN_Dec\t: '+str(np.min(healtbl['dec']))+'\n')
	f.close()
	os.system('chmod 777 '+save_path+'/*')
	#------------------------------------------------------------#
	#	6.	MAKE TARGET LIST FOR EACH OBSERVATORIES
	#------------------------------------------------------------#
	print('2-6.\tMAKE TARGET LIST FOR EACH OBSERVATORIES (.txt)')
	rtstbl		= Table()
	rtstbl['name'], rtstbl['ra'], rtstbl['dec'], rtstbl['sort'], rtstbl['rank'], rtstbl['dist']	= cantbl['iden_name'], cantbl['ra'], cantbl['dec'], cantbl['sort'], cantbl['rank'], cantbl['dist']
	rtstbl.write(save_path+'/'+eventname+'-rts_candi.txt', format='ascii', overwrite=True)
	obspath		= '/home/gw/Research/observatory.txt'
	catpath		= save_path+'/'+eventname+'-rts_candi.txt'
	if params['AlertType'] == 'Preliminary':
		deldays	= 0
	elif params['AlertType'] == 'Initial':
		deldays = 1
	elif params['AlertType'] == 'Update':
		deldays = 2
	else:
		deldays = 3
	# for obs in ['SAO']:
	for obs in ['SAO', 'McD', 'MAO', 'LOAO', 'SSO', 'UKIRT']:#, 'Spain', 'NMex', 'Calif']:
		print(obs)
		param_rts	= dict(	observatory=obs,
							headname=eventname,
							save_path=save_path,
							obspath=obspath,
							catpath=catpath,
							start=yyyy+'/'+mm+'/'+dd,
							end=yyyy+'/'+mm+'/'+str(int(dd)+deldays),
							altlimit=30.,
							moonseperation=40.,
							sunlimit='-18',
							numlimit=300)
		tool.rtsmaker(**param_rts)
	#------------------------------------------------------------#
	#	6.5 EXCEPTION FOR SQUEAN (REMOVE ':')
	#------------------------------------------------------------#
	'''
	rtsfile	= save_path+'/'+eventname+'-'+yyyy+mm+dd+'-rts_vis-McD.txt'
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
		print('{0}\n4-0.\tSEND MAIL TO CO-WORKERS{1}'.format('-'*60, len(matchtbl), '-'*60))
		reciver	= ascii.read('/home/gw/Research/alert-mail-reciver.txt')
		subject	= '[GW alert-{}] {}'.format(progenitor, eventname)
		# contents= 'CONTENTS'
		import codecs
		contents= codecs.open(save_path+'/'+eventname+'-summary.txt', 'rb', 'utf-8')

		fromID	= 'ceouobs@gmail.com'
		fromPW	= 'ceou@snu'
		toIDs	= ''
		for address in reciver['address']: toIDs += address+','
		toIDs	= toIDs[:-1]
		#toIDs	= "gregorypaek94@gmail.com"
		#toIDs	= 'joonhoyo305@gmail.com'
		ccIDs	= 'gundam_psh@naver.com'
		import glob
		path	= glob.glob(save_path+'/'+eventname+'-*.txt')
		#path.append(save_healfix+yymmdd+'/'+yymmdd+'-'+hhmmss+'-skymap.png')
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
		pnglist		= glob.glob(save_path+'/'+eventname+'*.png')
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
	# os.system('cat '+save_path+'/'+eventname+'-summary.txt')
	os.system('chmod 777 '+save_path)
	os.system('chmod 777 '+save_path+'/*')
	deltime		= time.time() - starttime
	print('DONE ({0} sec)'.format(round(deltime, 1)))
	print('MAIN PART FOR GW ALERT IS STANDBY.\n{0}\n\n\n'.format('='*60))
#------------------------------------------------------------#
#   Listen for GCNs until the program is interrupted
#   (killed or interrupted with control-C).
gcn.listen(handler=process_gcn)
#------------------------------------------------------------#