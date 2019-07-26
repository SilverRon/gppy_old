#   GW ALERT FOR ER 14 RUN IN QSO SERVER
#   BASED ON https://emfollow.docs.ligo.org/userguide/tutorial/index.html
#   2019.03.03 MADE		BY Gregory S.H. Paek
#	2019.03.28 UPDATED	BY Gregory S.H. Paek
#	2019.04.18 UPDATED	BY Gregory S.H. Paek
#	2019.05.02 UPDATED	BY Gregory S.H. Paek (ver. 2.0)
#	2019.05.10 UPDATED	BY Gregory S.H. Paek
#============================================================#
#	MODULE
#------------------------------------------------------------#
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
#------------------------------------------------------------#
starttime	= time.time()
catpath		= '/home/sonic/Research/cat/GLADE2.3/GLADE_2.3+2MASS_PSC+identi_name.dat'
#------------------------------------------------------------#
#	TEST BAYESTAR
#------------------------------------------------------------#
os.system('ls *.fits *fits.gz')
healpixfits	= input('BAYSTAR PATH\t: ')
# confidence	= float(input('CONFIDENCE (0.5,0.9)\t: '))
# role, eventname	= input('ROLE (test,observation)\t: '), input('EVENTNAME (test)\t: ')
# if confidence	== '': confidence = 0.9
# if role			== '': role = 'test'
# if eventname	== '': eventname = 'test'
confidence, role, eventname = 0.9, 'observation', 'S190425z_Update'

#------------------------------------------------------------#
save_healfix	= './'
save_path		= save_healfix+eventname
os.system('mkdir '+save_path)
yymmdd, hhmmss	= tool.timename()
yy	= yymmdd[0:2]
yyyy= '20'+yy
mm	= yymmdd[2:4]
dd	= yymmdd[4:6]
#------------------------------------------------------------#
#	1. bayestar.fits -> RA/Dec/P_2D TABLE
#------------------------------------------------------------#
comment     = '1.\tReading healpix file'; print(comment)
param_heal	= dict(	healpixfits=healpixfits,
					confidence=confidence,
					eventname=eventname,
					save_path=save_path,
					hdr=True,
					view=True)
healtbl, hdr= gw.heal2table(**param_heal)
healtbl.write(save_path+'/'+eventname+'_healpix.dat', format='ascii', overwrite=True)
gwdist      = hdr['DISTMEAN']
gwdiststd   = hdr['DISTSTD']
npix        = hdr['NAXIS2']         # i.e. 50331648
#------------------------------------------------------------#
#	2.	GLADE catalog -> LOAD AS ascii
#------------------------------------------------------------#
comment     = '2.\tReading GLADE catalog'; print(comment)
try:
	print(len(cattbl))
except:
	cattbl		= ascii.read(catpath)
#   CALCULATE MEAN DISTANCE BETW. PIXELS
alldeg		= (4*np.pi)*((180./np.pi)**2)         	#	[deg^2]
maxerr	    = 2.*np.sqrt(alldeg / npix)*3600.		#	["/pix]
#------------------------------------------------------------#
#	3.	MATCHING GLADE & HEALPIX
#------------------------------------------------------------#
comment		= '3.\tMATCHING START'; print(comment)
#	SKY COORD.
coo_intbl   = SkyCoord(healtbl['ra'], healtbl['dec'], unit=(u.deg, u.deg))
coo_reftbl  = SkyCoord(cattbl['ra'], cattbl['dec'], unit=(u.deg, u.deg))
#   INDEX FOR REF.TABLE
indx_match, d2d, d3d  = coo_intbl.match_to_catalog_sky(coo_reftbl, nthneighbor=1)
comment		= '\nMATHCING END\n'; print(comment)
#	REMOVE MATCHING & null & dist.cut
dummytbl		= cattbl[indx_match]
dummytbl['sep']	= d2d
dummytbl['P_2D']= healtbl['P_2D']
n=1
matchtbl0		= dummytbl[	(dummytbl['sep']*3600. < maxerr) &
							(dummytbl['dist'] > gwdist-n*gwdiststd) &
							(dummytbl['dist'] < gwdist+n*gwdiststd)]
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
#	4.	SEPERATE TABLE
#------------------------------------------------------------#
comment			= '4.\tSEPERATE TABLE'; print(comment)
#	K TABLE / NO K BUT B TABLE / NO K AND B TABLE (ONLY DIST)
cantblK			= matchtbl[	(matchtbl['K']!=-99)]
cantblB			= matchtbl[	(matchtbl['K']==-99) &
							(matchtbl['B']!=-99)]
cantbldist		= matchtbl[	(matchtbl['K']==-99) &
							(matchtbl['B']==-99)]
#------------------------------------------------------------#
#	5.	SORTING: P_2D & Luminosity & Distance
#------------------------------------------------------------#
comment			= '5.\tSORTING: P_2D & Luminosity & Distance'; print(comment)
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
		print('K-band info. is available.')
		cantblK['score']	= gw.prob(**param_K)
		cantblK				= cantblK[np.argsort(-1*cantblK['score'])]
		cantblK['sort']		= 'K'
		tbllist.append(cantblK)
	if	len(cantblB) > 0:
		print('B-band info. is available.')
		cantblB['score']	= gw.prob(**param_B)
		cantblB				= cantblB[np.argsort(-1*cantblB['score'])]
		cantblB['score']	= cantblB['score']*( np.min(cantblK['score'])*0.99/np.max(cantblB['score']) )
		cantblB['sort']		= 'B'
		tbllist.append(cantblB)
	if	len(cantbldist) > 0:
		print('Distance info. is available.')
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
	print('\nFAIL TO SORTING PROCESS. MAKE TABLE WITHOUT CONSTRAINTS.\n')
	cantbl		= Table(	[matchtbl['name'], matchtbl['ra'], matchtbl['dec']], matchtbl['P_2D']/np.sum(matchtbl['P_2D']),
							names=['name', 'ra', 'dec', 'score'])
	cantbl['sort']= 'none'

cantbl['rank']		= np.arange(len(cantbl))+1
cantbl.write(save_path+'/'+eventname+'-all_candi.txt', format='ascii', overwrite=True)
targettbl			= Table(	[cantbl['name'], cantbl['ra'], cantbl['dec'], cantbl['sort']],
								names=['name', 'ra', 'dec', 'sort'])
#------------------------------------------------------------#
#	PLOT CUMULATIVE SCORE
#------------------------------------------------------------#
plt.close('all')
plt.rcParams.update({'font.size': 20})
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
#------------------------------------------------------------#
#	MAKE REGION FILE
#------------------------------------------------------------#
param_regmaker	= dict(	filename	= save_path+'/'+eventname,
						intbl		= targettbl,
						name		= 'name',
						racol		= 'ra',
						deccol		= 'dec')
tool.ds9regmaker(**param_regmaker)
#------------------------------------------------------------#
#	WRITE SUMMARY
#------------------------------------------------------------#
f	= open(save_path+'/'+eventname+'-summary.txt', 'w')
skydeg_cut	= (alldeg/npix)*len(healtbl)
skydeg_cut_N= (alldeg/npix)*len(healtbl[healtbl['dec']>0])
skydeg_cut_S= (alldeg/npix)*len(healtbl[healtbl['dec']<0])
mag, magerr	= gw.kilonova_mag(gwdist, gwdiststd)
f.write('EXPECTED KILONOVA APP.MAG IN i-BAND [AB] : '+str(round(mag, 2))+'+/-'+str(round(magerr, 2))+'\n')
f.write(str(confidence*100)+'% confidence area\t: '+str(skydeg_cut)+'\t[deg^2]\n')
f.write(str(confidence*100)+'% confidence area for Northern sphere\t: '+str(skydeg_cut_N)+'\t[deg^2]\n')
f.write(str(confidence*100)+'% confidence area for Southern sphere\t: '+str(skydeg_cut_S)+'\t[deg^2]\n')
f.write('Number of galaxy candidates\t: '+str(len(matchtbl))+'\n')
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
comment		= '6.\tMAKE TARGET LIST FOR EACH OBSERVATORIES'; print(comment)
rtstbl		= Table()
rtstbl['name'], rtstbl['ra'], rtstbl['dec'], rtstbl['sort'], rtstbl['rank'], rtstbl['dist']	= cantbl['iden_name'], cantbl['ra'], cantbl['dec'], cantbl['sort'], cantbl['rank'], cantbl['dist']
rtstbl.write(save_path+'/'+eventname+'-rts_candi.txt', format='ascii', overwrite=True)

'''obspath		= '/home/gw/Research/observatory.txt'
catpath		= save_path+'/'+eventname+'-rts_candi.txt'
for obs in ['SAO', 'McD', 'MAO', 'LOAO', 'SSO', 'UKIRT']:#, 'Spain', 'NMex', 'Calif']:
	print(obs)
	param_rts	= dict(	observatory=obs,
						headname=eventname,
						save_path=save_path,
						obspath=obspath,
						catpath=catpath,
						start=yyyy+'/'+mm+'/'+dd,
						end=yyyy+'/'+mm+'/'+str(int(dd)+1),
						altlimit=30.,
						moonseperation=40.,
						sunlimit='-18',
						numlimit=300)
	tool.rtsmaker(**param_rts)'''
#------------------------------------------------------------#
#	7.	SEND MAIL
#------------------------------------------------------------#
'''
if (role != 'test') & (float(params['Terrestrial']) < 0.5):
	comment	= '7.\tSEND TARGET LIST TO CO-WORKERS'; print(comment)
	reciver	= ascii.read('/home/gw/Research/alert-mail-reciver.txt')
	subject	= '[GW alert-'+str(root.attrib['role'])+'] '+eventname+' GW EM follow-up observation'
	#contents= 'Dear co-workers,\nThis is test mail.\nminor detail will be updated.\n- Gregory S.H. Paek [+82 10-2285-8786]\n'
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
os.system('cat '+save_path+'/'+eventname+'-summary.txt')
os.system('chmod 777 '+save_path)
os.system('chmod 777 '+save_path+'/*')
deltime		= time.time() - starttime
print('All PROCESS IS DONE.\t('+str(round(deltime, 1))+' sec)')
print('MAIN PART FOR GW ALERT IS STANDBY.')
'''
