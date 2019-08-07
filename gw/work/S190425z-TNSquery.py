#============================================================
#	TNS QUERY CODE FOR S190425z WORK
#	2017.09.14	CREATED BY Nikola Knezevic
#	2019.08.06	MODIFIED BY Gregory S.H. Paek
#	2019.08.07	MODIFIED BY Gregory S.H. Paek
#============================================================
import os, glob
import numpy as np
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.wcs import WCS
import requests
import json
from collections import OrderedDict
from astropy.table import Table, vstack, Column, MaskedColumn
from imsng import phot
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
#	FUNCTION
#============================================================#
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
imtbl = ascii.read('tns.dat')
indexlist = [] 
for i in range(len(imtbl)): 
	trs = imtbl['transients'][i] 
	if '2019' in trs: 			#	SELECT RECENT ONES
		print(trs) 
		indexlist.append(i) 
seltbl = imtbl[indexlist]		#	SELECTED TABLE

for i in range(len(seltbl)):
	inim = seltbl['image'][i]
	w = WCS(inim)
	trs = seltbl['transients'][i]
	tnames, txs, tys = [], [], []
	for transient in trs.split(','):
		if '2019' in transient:
			qname = transient.split('_')[1]
			onetrtbl = query_transient_routine(qname)
			tra, tdec = onetrtbl['radeg'][0], onetrtbl['decdeg'][0]		#	FIRST ROW
			tx, ty = w.wcs_world2pix(tra, tdec, 0)
			tdateobs = onetrtbl['discoverydate'][0].replace(' ', 'T')
			ttime = Time(tdateobs, format='isot')
			tjd, tmjd = ttime.jd, ttime.mjd
			deljd = seltbl['jd'][i]-tjd
			if deljd > -3:
				tnames.append('{}_{}d_{}mag'.format(transient, round(deljd, 3), round(onetrtbl['discoverymag'][0], 3)))
				txs.append(tx)
				tys.append(ty)
			else:
				pass
		else:
			pass
		phot.plotshow(inim, tnames, txs, tys, outname=inim[:-5]+'_TNS.png')