#!/home/sonic/anaconda3/envs/astroconda/lib python3.7
# -*- coding: utf-8 -*-
#============================================================
#	TNS QUERY CODE
#	2017.09.14	CREATED BY Nikola Knezevic
#	2019.08.06	MODIFIED BY Gregory S.H. Paek
#============================================================
import os, glob
import numpy as np
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
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
	radius = np.asscalar(obstbl[obstbl['obs'] == obs]['fov'])   # ['], [arcmin]
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
def query_transient_routine(qname, url_tns_sandbox_api="https://sandbox-tns.weizmann.ac.il/api/get", photometry='1', spectra='1'):
	get_obj = [("objname", qname), ("photometry", photometry), ("spectra", spectra)]
	response = get(url_tns_sandbox_api, get_obj)
	if None not in response:
		cols = []
		json_data=format_to_json(response.text)
		json_dict_targ = json.loads(json_data)
		for key in ['objname', 'redshift', 'hostname', 'host_redshift', 'discoverydate', 'discoverymag', 'discmagfilter', 'name_prefix', 'type', 'discoverer', 'internal_name', 'ra', 'radeg', 'dec', 'decdeg', 'object_type']:
			mask = False
			try:
				if key == 'discmagfilter':
					val = json_dict_targ['data']['reply'][key]['name']
				elif key == 'object_type':
					val = json_dict_targ['data']['reply'][key]['name']
					val = val.replace(' ','')
				else:
					val = json_dict_targ['data']['reply'][key]
			except:
				val = -99
				mask = True
			column = MaskedColumn([val], name=key, mask=[mask])
			cols.append(column)
		intrtbl = Table( cols )
		return intrtbl
	else:
		print (response[1])
		return response[1]
#============================================================#
#	MAIN BODY
#============================================================#
path_and_file='/home/sonic/Research/table/obs.txt'
obstbl = ascii.read(path_and_file)
imlist = glob.glob('*.fits')

n = 0
imtblist = []
for inim in imlist:
	n += 1
	print('PROCESS [{}/{}]\t: {}'.format(n, len(imlist), inim))
	oneimtbl = search_transient_routine(inim)
	transients = np.asscalar(oneimtbl['transients'])
	imtblist.append(oneimtbl)
	trtblist = []
	if transients != 'None':
		for transient in transients.split(','):
			qname = transient.split('_')[1]		#	QUERY NAME
			onetrtbl = query_transient_routine(qname)
			trtblist.append(onetrtbl)
		try:
			trtbl = vstack(trtblist)
			trtbl.write('{0}_TNS.dat'.format(inim[:-5]), format='ascii', overwrite=True)
		except:
			print('FAIL')
			pass
	else:
		pass

imtbl = vstack(imtblist)
imtbl.write('tns.dat', format='ascii', overwrite=True)


"""
# get ascii file (according to the "asciifile" name obtained from the get obj reply)
ascii_file_url="https://sandbox-tns.weizmann.ac.il/system/files/uploaded/"\
			   "Padova-Asiago/tns_2017A_2457777.69_Ekar_AFOSC_Padova-Asiago.txt"
get_file(ascii_file_url)

# get fits file (according to the "fitsfile" name obtained from the get obj reply)
fits_file_url="https://sandbox-tns.weizmann.ac.il/system/files/uploaded/"\
			  "Padova-Asiago/tns_2017A_2457777.69_Ekar_AFOSC_Padova-Asiago.fits"
get_file(fits_file_url)

"""




