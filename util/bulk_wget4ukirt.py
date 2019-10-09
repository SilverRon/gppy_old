#!/home/sonic/anaconda3/envs/astroconda/bin/python
# -*- coding: utf-8 -*-
#============================================================
#	wget UKIRT DATA
#	2019.09.26	CREATED BY Gregory S.H. Paek
#============================================================
import os, glob
import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table, vstack, Column, MaskedColumn
#============================================================
#	FUNCTION
#------------------------------------------------------------
def wgetcom(user, password, path_server, path_save, filename):
	com = 'wget --user={0} --password={1} -r {2}/{3} -P {4}'.format(user, password, path_server, filename, path_save)
	print(com)
	os.system(com)
	# com = 'wget --user=u19aeap001 --password=090185 -r http://apm3.ast.cam.ac.uk/~mike/wfcam/u19aeap001/w20190913_00385_sf_st.fit -P /data1/UKIRT'
#------------------------------------------------------------
path_key = '/home/sonic/Research/table/keys.dat'
path_server = 'http://apm3.ast.cam.ac.uk/~mike/wfcam/u19aeap001'
path_save = '/data1/UKIRT'
proj_name = 'ukirt19a'
#------------------------------------------------------------
os.system('ls /data1/UKIRT/*')
path_table = input('WHICH DATA TABLE :\t')
#------------------------------------------------------------
keys = ascii.read(path_key)
proj_idx = np.where(keys['name']==proj_name)
proj_id, proj_pw = np.asscalar(keys['key'][proj_idx]), np.asscalar(keys['pw'][proj_idx])
proj_pw = str(proj_pw)
if len(proj_pw)<=6:
	proj_pw = '0'+proj_pw
#------------------------------------------------------------
#	DOWNLOAD
#------------------------------------------------------------
intbl = ascii.read(path_table)
for i, filename in enumerate(intbl['image']):
	print('[{}/{}]'.format(i+1, len(intbl)))
	param_wgetcom = dict(	user = proj_id,
							password = proj_pw,
							path_server = path_server,
							path_save = path_save,
							filename = filename)
	wgetcom(**param_wgetcom)

path_part = [path_save+'/', path_server.split('//')[1]+'/', '*']
path_files = ''.join(path_part)
os.system('mv {} {}'.format(path_files, path_save))
os.system('rm -rf {}'.format(path_files[:-2]))