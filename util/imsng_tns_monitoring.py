#!/home/sonic/anaconda3/envs/astroconda/lib python3.7
# -*- coding: utf-8 -*-
#============================================================
#	IMSNG MONITORING NEW TRANSIENT WITH TNS QUERY CODE
#	2017.09.14	CREATED BY Nikola Knezevic
#	2019.08.06	MODIFIED BY Gregory S.H. Paek
#	2019.08.07	MODIFIED BY Gregory S.H. Paek
#============================================================
import os, glob
import numpy as np
import time, datetime
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import requests
import json
from collections import OrderedDict
from astropy.table import Table, vstack, Column, MaskedColumn
############################# PARAMETERS #############################
# API key for Bot                                                    #
#	GET KEY FROM PERSONAL TABLE
path_keys = '/home/sonic/Research/table/keys.dat'
keytbl = ascii.read(path_keys)
api_key = np.asscalar(keytbl['key'][keytbl['name']=='TNS'])          #
# list that represents json file for search obj                      #
search_obj=[("ra",""), ("dec",""), ("radius",""), ("units",""),      #
			("objname",""), ("internal_name","")]                    #
# list that represents json file for get obj                         #
get_obj=[("objname",""), ("photometry","0"), ("spectra","1")]        #
######################################################################

#############################    URL-s   #############################
# url of TNS and TNS-sandbox api                                     #
url_tns_api="https://wis-tns.weizmann.ac.il/api/get"                 #
url_tns_sandbox_api="https://sandbox-tns.weizmann.ac.il/api/get"     #
######################################################################

############################# DIRECTORIES ############################
# current working directory                                          #
cwd=os.getcwd()                                                      #
# directory for downloaded files                                     #
download_dir=os.path.join(cwd,'downloaded_files')                    #
######################################################################

########################## API FUNCTIONS #############################
# function for changing data to json format                          #
def format_to_json(source):                                          #
	# change data to json format and return                          #
	parsed=json.loads(source,object_pairs_hook=OrderedDict)          #
	result=json.dumps(parsed,indent=4)                               #
	return result                                                    #
#--------------------------------------------------------------------#
# function for search obj                                            #
def search(url,json_list):                                           #
	try:                                                               #
		# url for search obj                                             #
		search_url=url+'/search'                                         #
		# change json_list to json format                                #
		json_file=OrderedDict(json_list)                                 #
		# construct the list of (key,value) pairs                        #
		search_data=[('api_key',(None, api_key)),                        #
					('data',(None,json.dumps(json_file)))]              #
		# search obj using request module                                #
		response=requests.post(search_url, files=search_data)            #
		# return response                                                #
		return response                                                  #
	except Exception as e:                                             #
		return [None,'Error message : \n'+str(e)]                        #
#--------------------------------------------------------------------#
# function for get obj                                               #
def get(url,json_list):                                              #
	try:                                                               #
		# url for get obj                                                #
		get_url=url+'/object'                                            #
		# change json_list to json format                                #
		json_file=OrderedDict(json_list)                                 #
		# construct the list of (key,value) pairs                        #
		get_data=[('api_key',(None, api_key)),                           #
					('data',(None,json.dumps(json_file)))]              #
		# get obj using request module                                   #
		response=requests.post(get_url, files=get_data)                  #
		# return response                                                #
		return response                                                  #
	except Exception as e:                                             #
		return [None,'Error message : \n'+str(e)]                        #
#--------------------------------------------------------------------#
# function for downloading file                                      #
def get_file(url):                                                   #
	try:                                                             #
		# take filename                                                  #
		filename=os.path.basename(url)                                   #
		# downloading file using request module                          #
		response=requests.post(url, files=[('api_key',(None, api_key))], #
							stream=True)                              #
		# saving file                                                    #
		path=os.path.join(download_dir,filename)                         #
		if response.status_code == 200:                                  #
			with open(path, 'wb') as f:                                  #
				for chunk in response:                                   #
					f.write(chunk)                                       #
			print ('File : '+filename+' is successfully downloaded.')    #
		else:                                                            #
			print ('File : '+filename+' was not downloaded.')            #
			print ('Please check what went wrong.')                      #
	except Exception as e:                                               #
		print ('Error message : \n'+str(e))                              #
######################################################################

#============================================================#
#	FUNCTION FOR ROUTINE
#============================================================#
def search_transient_routine(inim, url_tns_sandbox_api="https://sandbox-tns.weizmann.ac.il/api/get"):
	obs = inim.split('-')[1]
	# obs = 'KMTNET'
	hdr = fits.getheader(inim)
	dateobs = hdr['DATE-OBS']
	t = Time(dateobs, format='isot')
	jd = t.jd
	mjd = t.mjd
	radeg, dedeg = hdr['CRVAL1'], hdr['CRVAL2']
	c = SkyCoord(radeg, dedeg, unit='deg')
	radec = c.to_string('hmsdms', sep=':')
	ra, dec = radec.split(' ')[0], radec.split(' ')[1]
	radius = np.asscalar(obstbl[obstbl['obs'] == obs]['fov'])/2   # ['], [arcmin]
	units = 'arcmin'
	#------------------------------------------------------------
	#	SEARCH OBJECT
	#------------------------------------------------------------
	search_obj=[("ra",ra), ("dec",dec), ("radius",radius), ("units",units),
				("objname",""), ("internal_name","")]                    
	response=search(url_tns_sandbox_api,search_obj)
	tnames = []
	if None not in response:
		json_data =format_to_json(response.text)
		json_dict = json.loads(json_data)
		if len(json_dict['data']['reply']) != 0:
			transients = ''
			for i in range(len(json_dict['data']['reply'])):
				tname = json_dict['data']['reply'][i]['prefix']+'_'+json_dict['data']['reply'][i]['objname']
				tnames.append(tname)
				transients = transients+tname+','
			transients = transients[:-1]
		else:
			transients = 'None'
	else:
		transients = 'None'
		print(response[1])
	#------------------------------------------------------------
	inimtbl = Table(	[[inim], [round(radeg, 5)], [round(dedeg, 3)], [ra], [dec], [dateobs], [round(jd, 5)], [round(mjd, 5)], [transients]],
					names=('image', 'radeg', 'dedeg', 'ra', 'dec', 'dateobs', 'jd', 'mjd', 'transients'))
	return inimtbl
#------------------------------------------------------------
def query_transient_routine_simple(qname, field, c, url_tns_sandbox_api="https://sandbox-tns.weizmann.ac.il/api/get", photometry='1', spectra='1'):
	get_obj = [("objname", qname), ("photometry", photometry), ("spectra", spectra)]
	response = get(url_tns_sandbox_api, get_obj)
	#------------------------------------------------------------
	if None not in response:
		rows = [field]
		cols = ['field']
		json_data=format_to_json(response.text)
		json_dict_targ = json.loads(json_data)
		#------------------------------------------------------------
		#	EXTRACT INFORMATION
		#------------------------------------------------------------		
		tc = SkyCoord(json_dict_targ['data']['reply']['ra'],
		json_dict_targ['data']['reply']['dec'], unit=(u.hourangle, u.deg))
		for key in ['objname', 'discoverydate', 'discoverymag', 'discmagfilter', 'internal_name', 'ra', 'radeg', 'dec', 'decdeg', 'object_type']:
			if key == 'discmagfilter':
				val = json_dict_targ['data']['reply'][key]['name']
			elif key == 'object_type':
				val = json_dict_targ['data']['reply'][key]['name']
			else:
				val = json_dict_targ['data']['reply'][key]
			rows.append(val) 
			cols.append(key)
		try:
			rows.append(round(c.separation(tc).arcminute, 3))
		except:
			rows.append(-99)
		cols.append('sep_arcmin')
		#------------------------------------------------------------
		rows = tuple(rows)
		cols = tuple(cols)
		return rows
	else:
		print (response[1])
		return response[1]
#------------------------------------------------------------
def getCurrentStrTime():
   return time.strftime("%Y%m%d-%H%M%S")
#============================================================#
#	MAIN BODY
#============================================================#
for i in range(206265):
	starttime = time.time()
	outname = 'IMSNG-TNS-{}.dat'.format(getCurrentStrTime())
	#------------------------------------------------------------
	path_obs = '/home/sonic/Research/table/obs.txt'
	path_input = '/home/sonic/Research/table/imsng-alltarget.dat'
	path_save = '/data1/IMSNG'
	refcats = glob.glob('/data1/IMSNG/*.dat')
	refcats.sort()
	path_ref = refcats[-1]
	#------------------------------------------------------------
	obstbl = ascii.read(path_obs)
	reftbl = ascii.read(path_ref)
	intbl = ascii.read(path_input)
	#------------------------------------------------------------
	radius = 15
	units = 'arcmin'
	#------------------------------------------------------------
	cols = ('field',
			'objname',
			'discoverydate',
			'discoverymag',
			'discmagfilter',
			'internal_name',
			'ra',
			'radeg',
			'dec',
			'decdeg',
			'object_type',
			'sep_arcmin')
	#------------------------------------------------------------
	tblist = []
	rowlist = []
	for i in range(len(intbl)):
		obj = intbl['obj'][i]
		ra, dec = intbl['ra'][i], intbl['dec'][i]
		c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

		print('PROCESS\t{}\t[{}/{}]'.format(obj, i+1, len(intbl)))
		search_obj = [("ra",ra), ("dec",dec), ("radius",radius), ("units",units),
					("objname",""), ("internal_name","")]                    
		response = search(url_tns_sandbox_api,search_obj)
		tnames = []
		if None not in response:
			json_data = format_to_json(response.text)
			json_dict = json.loads(json_data)
			if len(json_dict['data']['reply']) != 0:
				transients = ''
				for i in range(len(json_dict['data']['reply'])):
					tname = json_dict['data']['reply'][i]['objname']
					tnames.append(tname)
					transients = transients+tname+','
				transients = transients[:-1]
			else:
				transients = 'None'
		else:
			transients = 'None'
			print(response[1])
		for qname in tnames:
			# onetbl = query_transient_routine_simple(qname)
			# onetbl['imsng'] = intbl['obj'][i]
			# tblist.append(onetbl)
			onerows = query_transient_routine_simple(qname, field=obj, c=c)
			rowlist.append(onerows)
	comtbl = Table(rows=rowlist, names=cols)
	#------------------------------------------------------------
	#	CHECK NEW TRANSIENTs
	#------------------------------------------------------------
	newlist = []
	for i, objname in enumerate(comtbl['objname']):
		if objname not in reftbl['objname']:
			newlist.append(comtbl[i])
	if len(newlist) == 0:
		print('THERE IS NO NEW ONE.')
		pass
	else:
		print('NEW TRANSIDENT WAS REPORTED. SEND E-MAIL TO STAFFs...')
		#------------------------------------------------------------
		#	SAVE TABLEs
		#------------------------------------------------------------
		comtbl.write(path_save+'/'+outname, format='ascii', overwrite=True)
		newtbl = vstack(newlist)
		newtbl.write(path_save+'/NEW-'+outname, format='ascii', overwrite=True)
		#------------------------------------------------------------
		#	MAIL SETTING
		#------------------------------------------------------------
		reciver = ascii.read('/home/sonic/Research/table/imsng-mail-reciver.txt')
		subject	= '[IMSNG] {} NEW TRANSIENTs'.format(getCurrentStrTime())
		# contents= 'CONTENTS'
		import codecs
		contents= codecs.open(path_save+'/NEW-'+outname, 'rb', 'utf-8')
		fromID	= 'ceouobs@gmail.com'
		fromPW	= 'ceou@snu'
		toIDs	= ''
		for address in reciver['address']: toIDs += address+','
		toIDs	= toIDs[:-1]
		#toIDs	= "gregorypaek94@gmail.com"
		# ccIDs	= 'gundam_psh@naver.com'
		import glob
		# path	= glob.glob(save_path+'/'+eventname+'-*.txt')
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
		msg.attach(MIMEText(contents.read()))
		'''
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
		'''
		#------------------------------------------------------------
		#	SEND MAIL
		#------------------------------------------------------------
		#	ACCESS TO GMAIL
		smtp_gmail = smtplib.SMTP_SSL('smtp.gmail.com', 465)
		smtp_gmail.login(fromID, fromPW)
		# smtp_gmail.sendmail(msg["From"], msg["To"].split(",") + msg["Cc"].split(","), msg.as_string())
		smtp_gmail.sendmail(msg["From"], msg["To"].split(","), msg.as_string())
		smtp_gmail.quit()
		# comment	= 'Send '+str(path)+'\nFrom\t'+fromID+'\nTo'; print(comment); print(toIDs)
		print('SENDING E-MAIL COMPLETE.')

	deltime	= time.time() - starttime
	print('\nDone.\t\t[{0} sec]'.format(round(deltime, 1)))
	time.sleep(3600)