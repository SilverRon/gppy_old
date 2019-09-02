#	USEFUL FUNC.(TOOL) IN IMSNG MODULE
#	2019.03.03	CREATED BY Gregory S.H. Paek
#	2019.08.29	UPDATED BY Gregory S.H. Paek
#============================================================
def timename():
	'''
	CONVERT 'TIME' TO YYMMDD, HHMMSS FORM.
	INPUT	:	NONE
	OUTPUT	:	STRIG FORM OF 'YYMMDD', 'HHMMSS'
	'''
	import numpy as np
	import time
	now				= time.gmtime(time.time())
	y, m, d			= now.tm_year, now.tm_mon, now.tm_mday
	ho, mi, se		= now.tm_hour, now.tm_min, now.tm_sec
	yy				= str(y)[2:]
	if len(str(m)) < 2:
		mm			= '0'+str(m)
	else:
		mm			= str(m)
	if len(str(d)) < 2:
		dd			= '0'+str(d)
	else:
		dd			= str(d)
	if len(str(ho)) < 2:
		hour		= '0'+str(ho)
	else:
		hour		= str(ho)
	if len(str(mi)) < 2:
		mini		= '0'+str(mi)
	else:
		mini		= str(mi)
	if len(str(se)) < 2:
		sec			= '0'+str(se)
	else:
		sec			= str(se)
	yymmdd			= yy+mm+dd
	hhmmss			= hour+mini+sec
	return yymmdd, hhmmss
#------------------------------------------------------------
def detection(name, ra, dec, time, location):
	import numpy as np
	import os, glob, sys
	from astropy import units as u
	from astropy.time import Time
	from astropy.coordinates import SkyCoord, EarthLocation, AltAz
	from astropy.coordinates import get_sun, get_moon
	from astropy.io import ascii
	from astropy.table import Table, Column
	target      = SkyCoord(ra, dec, unit='deg') # defaults to ICRS frame
	site        = location
	del_midnight= np.linspace(-12, +12, 720) * u.hour
	time_night  = time+del_midnight
	frame_night = AltAz(obstime=time_night, location=site)

	targetaltaz_night   = target.transform_to(frame_night)
	sunaltaz_night      = get_sun(time_night).transform_to(frame_night)

	# indx_set            = np.where( sunaltaz_night.alt > -18 * u.deg )
	indx_rise           = np.where( sunaltaz_night.alt < -18 * u.deg )
	
	sunset              = del_midnight[np.min(indx_rise)]
	sunrise             = del_midnight[np.max(indx_rise)]
	
	del_midnight= np.linspace(sunset.value, sunrise.value, 100) * u.hour
	time_night  = time+del_midnight
	frame_night = AltAz(obstime=time_night, location=site)
	targetaltaz_night   = target.transform_to(frame_night)

	return targetaltaz_night
#------------------------------------------------------------
def sendmail(filename, subject, sendID, sendPW, reciver):
	'''
	Security reference
	https://cpuu.postype.com/post/23066
	Code reference
	https://kimdoky.github.io/python/2017/07/21/smtplib_email.html
	File attach
	https://brunch.co.kr/@jk-lab/31
	'''
	import smtplib
	from email.mime.text import MIMEText
	import codecs
	email_text = codecs.open(filename, 'rb', 'utf-8')
	msg = MIMEText(email_text.read())
	email_text.close()
	msg['Subject']	= subject
	msg['From']		= sendID
	smtp_gmail		= smtplib.SMTP_SSL('smtp.gmail.com', 465)
	smtp_gmail.login(sendID, sendPW)
	smtp_gmail.sendmail(sendID, reciver, msg.as_string())
	smtp_gmail.quit()
	comment	= 'Send '+filename+'\n'+'From '+sendID+' To '+reciver; print(comment)
#------------------------------------------------------------
def send_gmail(subject, contents, fromID, fromPW, toIDs, ccIDs=None, path=None):
	'''
	SEND GMAIL
	Security reference
	https://cpuu.postype.com/post/23066
	Code reference
	https://kimdoky.github.io/python/2017/07/21/smtplib_email.html
	File attach
	https://brunch.co.kr/@jk-lab/31
	'''
	import os
	import smtplib
	from email.mime.base import MIMEBase
	from email.mime.text import MIMEText
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
	msg.attach(MIMEText(contents, 'plain', 'utf-8'))
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
	


	#	ACCESS TO GMAIL & SEND MAIL
	smtp_gmail		= smtplib.SMTP_SSL('smtp.gmail.com', 465)
	smtp_gmail.login(fromID, fromPW)
	smtp_gmail.sendmail(msg["From"], msg["To"].split(",") + msg["Cc"].split(","), msg.as_string())
	smtp_gmail.quit()
	comment	= 'Send '+str(path)+'\nFrom\t'+fromID+'\nTo'; print(comment); print(toIDs)
#------------------------------------------------------------
def abs2app(mag, magerr, gwdist, gwdiststd):
	import numpy as np
	app		= 5*np.log10(gwdist)-5+mag
	apperr	= 5*gwdiststd/(gwdist*np.log(10))
	return app, apperr
#------------------------------------------------------------
def GW170817_like(gwdist, gwdiststd):
	import numpy as np
	m0		= 17.476	#	[AB] in i-band (t0+10h)
	m0err	= 0.018	
	dist0	= 38.4		#	[MPC]	Im et al. 2017
	dist0err= 8.9
	m		= m0+5.*np.log10(gwdist/dist0)
	merr	= np.sqrt( (m0err)**2 + ((5.*gwdiststd)/(gwdist*np.log(10)))**2 + ((5.*dist0err)/(dist0*np.log(10)))**2 )
	return m, merr
#------------------------------------------------------------
def func_linear(a, x, scaling=[0, 0]):
	xpt, ypt= scaling[0], scaling[1]
	ydel	= ypt - (-1*a*xpt)
	return -1*a*x + ydel
#------------------------------------------------------------
def calc_app(mag, magerr, gwdist0, gwdiststd0, gwdist1, gwdiststd1):
	import numpy as np
	app		= mag+5*np.log10(gwdist1/gwdist0)
	apperr	= np.sqrt( (magerr)**2 + ((5*gwdiststd1)/(np.log(5)*gwdist1))**2 + ((5*gwdiststd0)/(np.log(5)*gwdist0))**2 )
	return app, apperr
#------------------------------------------------------------
def ds9regmaker(filename, name, ra, dec):
	import os,sys
	import string
	from astropy.io import ascii 
	import numpy as np
	import math
	'''
	racol		= 'ALPHA_J2000'
	deccol		= 'DELTA_J2000'
	name		= 'NUMBER'
	intbl		= ascii.read(filename)
	'''
	radius	= """ 5" """
	color	= "green"
	f		= open(filename, 'w')

	head1	= "# Region file format: DS9 version 4.1\n"
	head2	= """global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n"""
	head3	= "fk5\n"

	f.write(head1)
	f.write(head2)
	f.write(head3)

	for n in range(len(ra)):
		body="circle("+str(ra[n])+","+str(dec[n])+","+radius+") # color="+color+" text={"+str(name[n])+"}\n"	
		f.write(body)
	f.close()
#------------------------------------------------------------
def rtsmaker(observatory, headname, save_path, obspath, catpath, start, end, altlimit=30., moonseperation=40., sunlimit='-18', numlimit=100):
	import pytz
	import jdcal
	import ephem
	#from numpy import *
	import numpy as np
	import os, sys
	import string
	import datetime
	import astropy.units as u
	from astropy.io import ascii
	import mskpy.observing as obs
	import astropy.coordinates as coord
	from astropy import units as u
	from astropy.coordinates import SkyCoord

	#------------------------------------------------------------#
	#	INPUT SAMPLE
	#------------------------------------------------------------#
	'''
	observatory		= 'SAO'
	save_path		= './'
	obspath			= "/home/gw/Research/observatory.txt"
	catpath			= 'MS181101ab_Preliminary-all_candidates.txt'
	start			= '2019/04/17'
	end				= '2019/04/19'
	#altitute limit and moon seperation, moon serperation is a little bit close (2~3 deg)
	numlimit		= 100
	altlimit		= 30.
	moonseperation	= 40.
	sunlimit		= '-18'
	'''
	#------------------------------------------------------------#
	#	OBSERVATORY INFO.
	#------------------------------------------------------------#
	obsinfo			= ascii.read(obspath)
	obsname     	= np.copy(obsinfo['name'])
	obsindex		= np.where(obsname == observatory)[0]
	obslat			= (np.copy(obsinfo['latitude(N+)'])[obsindex])[0]
	obslon			= (np.copy(obsinfo['longitude(E+)'])[obsindex])[0]
	obsalt			= (np.copy(obsinfo['altitude'])[obsindex])[0]
	obstz			= (np.copy(obsinfo['timezone'])[obsindex])[0]
	tz				= pytz.timezone(obstz)
	#------------------------------------------------------------#
	observ			= ephem.Observer()
	observ.lat		= str(obslat)
	observ.lon		= str(obslon)
	observ.elevation= obsalt
	observ.horizon	= sunlimit
	#------------------------------------------------------------#
	#objects from catalog file
	tdata			= ascii.read(catpath)
	objname			= tdata['name']
	ra				= tdata['ra']
	dec				= tdata['dec']
	prior			= tdata['sort']
	rank			= tdata['rank']			
	dist			= tdata['dist']

	RA				= coord.Angle(ra, unit = u.deg)
	Dec				= coord.Angle(dec, unit = u.deg)

	radd			= RA.value
	rad				= RA.hour
	decd			= Dec.value
	decdd			= Dec.degree

	#angular distance calculation
	def angsep(ra1deg, dec1deg, ra2deg, dec2deg) : 
		ra1rad		= ra1deg*np.pi/180
		dec1rad		= dec1deg*np.pi/180
		ra2rad		= ra2deg*np.pi/180
		dec2rad		= dec2deg*np.pi/180
		cos_a		= np.sin(dec1rad)*np.sin(dec2rad)+(np.cos(dec1rad)*np.cos(dec2rad)*np.cos(ra1rad-ra2rad))
		anglesep	= np.arccos(cos_a)*180/np.pi
		return anglesep

	#dates to calculate
	fmt				= '%Y/%m/%d'
	startdt			= datetime.datetime.strptime(start, fmt)
	enddt			= datetime.datetime.strptime(end, fmt)
	startmjd		= (jdcal.gcal2jd(startdt.year, startdt.month, startdt.day))[1]
	endmjd			= (jdcal.gcal2jd(enddt.year, enddt.month, enddt.day))[1]

	for i in range(int(endmjd-startmjd+1)):
		onedaymjd   = startmjd+i+1
		oneday      = jdcal.jd2gcal(2400000.5, onedaymjd)
		onedaydt    = datetime.datetime(oneday[0], oneday[1], oneday[2])
		dst         = tz.dst(onedaydt, is_dst=True)
		dst         = dst.seconds/3600

		onedaydt    = datetime.datetime(oneday[0], oneday[1], oneday[2], tzinfo=tz)
		onedayutc   = onedaydt.astimezone(pytz.utc)
		observ.date = onedayutc

		# Moon distance and information
		mcoord		= ephem.Moon()
		mcoord.compute(observ)
		minfo		= 'Moon ra, dec : '+str(mcoord.ra)+' '+str(mcoord.dec)+'\n'
		mphase		= ephem.Moon(observ.date)
		mphasestr	= 'Moon phase   : '+ "%.2f" % mphase.moon_phase +'\n'
		msep		= angsep(radd, decdd, np.degrees(mcoord.ra), np.degrees(mcoord.dec))

		#	SUNSET CALC.
		sunset      = observ.previous_setting(ephem.Sun())
		sunsettu    = ephem.Date.tuple(sunset)
		sunsetdt    = datetime.datetime(sunsettu[0],sunsettu[1],sunsettu[2],sunsettu[3],int(sunsettu[4]),tzinfo=pytz.utc)
		sunsetlocal = sunsetdt.astimezone(tz)
		sunsetstr   = sunlimit+' deg sunset : '+str(sunsetlocal.hour)+':'+str(sunsetlocal.minute)+'\n'
		sunsethour  = sunsetlocal.hour+sunsetlocal.minute/60.+sunsetlocal.second/3600.
		#	SUNRISE CALC.
		sunrise      = observ.next_rising(ephem.Sun())
		sunrisetu    = ephem.Date.tuple(sunrise)
		sunrisedt    = datetime.datetime(sunrisetu[0],sunrisetu[1],sunrisetu[2],sunrisetu[3],int(sunrisetu[4]),tzinfo=pytz.utc)
		sunriselocal = sunrisedt.astimezone(tz)
		sunrisestr   = sunlimit+' deg sunrise : '+str(sunriselocal.hour)+':'+str(sunriselocal.minute)+'\n'
		sunrisehour  = sunriselocal.hour+sunriselocal.minute/60.+sunriselocal.second/3600.
		#print (observatory)
		#print ('Local mid night in UTC : '+str(observ.date))
		#print (minfo,mphasestr,sunsetstr,sunrisestr)

		#	MAKE RESULT FILE
		stryear      = str(oneday[0])
		strmonth     = str(oneday[1])
		strday       = str(oneday[2]-1)
		if int(strmonth) < 10 : strmonth = '0'+strmonth
		if int(strday)   < 10 : strday   = '0'+strday

		f			= open(save_path+'/'+headname+'-'+stryear+strmonth+strday+"-rts_vis-"+observatory+".txt",'w')
		f.write('#\t'+str(observ.date)+' UTC & Day Time Saving +'+str(dst)+'\n')
		f.write('#\tObservatory\t= '+observatory+'\n')
		f.write('#\t'+sunsetstr)
		f.write('#\t'+sunrisestr)
		f.write('#\t'+minfo)
		f.write('#\t'+mphasestr)
		f.write('#\tMoon seperation = '+str(moonseperation)+'\n')
		f.write('#\tAltitude limit = '+str(altlimit)+'\n')
		f.write('#\tRank : the lower rank, the higher priority\n')
		f.write('#------------------------------------------------------- \n')
		f.write('name ra dec rise(LT) transit(LT) set(LT) moon_dist(deg) distance(Mpc) rank\n')

		numcount	= 0
		for n in range(len(rad)):
			#calculate rise transit set time with altitute limit
			param_rts	= dict(	ra=radd[n],
								dec=decdd[n],
								date=onedaydt,
								lon=obslon,
								lat=obslat,
								tz=obstz,
								limit=altlimit,
								precision=1440)
			rtscal		= obs.rts(**param_rts)
			rt			= rtscal[0]
			tt			= rtscal[1]
			st			= rtscal[2]
			if rtscal[0]== None:
				#print (objname[n],ra[n],dec[n], rtscal[0], rtscal[1], rtscal[2],"%.2f" % msep[n])
				pass
			elif sunrisehour < rtscal[0] < sunsethour and sunrisehour < rtscal[2] < sunsethour and sunrisehour < rtscal[1] < sunsethour: 
				#print (objname[n]+' It can be seen in daytime!')
				pass
			elif msep[n] < moonseperation or msep[n] > 360-moonseperation:
				#print (objname[n]+' too close to Moon < '+str(moonseperation)+' deg')
				pass
			else:
				if numcount < numlimit:
					c= SkyCoord(ra=ra[n]*u.degree, dec=dec[n]*u.degree, frame='icrs')
					c_ra= c.ra.hms
					c_dec= c.dec.dms
					
					nra='%02d:%02d:%.3f' %(c_ra[0], abs(c_ra[1]), abs(c_ra[2]))
					ndec='%02d:%02d:%.3f' %(c_dec[0], abs(c_dec[1]), abs(c_dec[2]))
					
					rtp	="%.2d" % int(rt)+':'+"%.2d" % int((rt-int(rt))*60)
					ttp	="%.2d" % int(tt)+':'+"%.2d" % int((tt-int(tt))*60)
					stp	="%.2d" % int(st)+':'+"%.2d" % int((st-int(st))*60)
					vis	='{:8s}	{:12s}	{:12s}	{:5s}	{:5s}	{:5s}	{:3s}		{:4s}	{:4s}'.format(objname[n],str(nra),str(ndec),rtp,ttp,stp,str(int(msep[n])),str(int(dist[n])),str(rank[n]))+'\n'
					f.write(vis)
					#print (objname[n],ra[n],dec[n], rtp,ttp,stp,"%.2f" % msep[n])
					numcount+= 1
				else:
					pass
				'''
				if numcount < numlimit:
					
					rtp="%.2d" % int(rt)+':'+"%.2d" % int((rt-int(rt))*60)
					ttp="%.2d" % int(tt)+':'+"%.2d" % int((tt-int(tt))*60)
					stp="%.2d" % int(st)+':'+"%.2d" % int((st-int(st))*60)
					vis='{:24s}	{:12s}	{:12s}	{:5s}	{:5s}	{:5s}	{:3s}	{:2s}'.format(objname[n],str(ra[n]),str(dec[n]),rtp,ttp,stp,str(int(msep[n])),str(prior[n]))+'\n'
					f.write(vis)
					#print (objname[n],ra[n],dec[n], rtp,ttp,stp,"%.2f" % msep[n])
					numcount+= 1
				else:
					pass
				'''
		f.close()

