#============================================================
#	PHOTOMETRY + SCORE FOR S190425z
#	2019.08.11	CREATED BY	Gregory S.H. Paek
#============================================================
import os, glob
import numpy as np
from astropy.io import ascii
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
from astropy.time import Time
from scipy.interpolate import interp1d
from astropy.table import Column
#============================================================
#	FUNCTION
#------------------------------------------------------------
def countscore(objlist, cantbl):
	score = 0.0
	for obj in objlist:
		indx = np.where(obj == cantbl['name'])
		score += np.asscalar(cantbl['score'][indx])
	return score
#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_phot = '/data1/S190425z/1.result/Initial/phot_ini_sources_GECKO.dat'
path_cand = '/data1/S190425z/info/Initial/S190425z_Initial-all_candi.txt'
path_save = '/data1/S190425z/1.result/Initial'
#------------------------------------------------------------
photbl = ascii.read(path_phot)
cantbl = ascii.read(path_cand)
obslist = ['KMTNET', 'LOAO', 'SAO', 'SQUEAN', 'LSGT', 'UKIRT']
#------------------------------------------------------------
#	SUMMARY
#------------------------------------------------------------
f = open(path_save+'/result4initial.dat', 'a')
f.write('#obs numb score set_numb set_score\n')
objslist = []
# scorelist = []
totscore = 0.0
totnumb = 0.0
set_totscore = 0.0
set_totnumb = 0.0
#------------------------------------------------------------
for obs in obslist:
	obstbl = photbl[photbl['obs']==obs]
	objlist = []
	# score = 0.0
	for src in obstbl['sources']:
		if src == 'None':
			pass
		elif ',' in src:	#	MULTIPLE SOURCES
			for sr in src.split(','):
				objlist.append(sr)
				# indx = np.where(sr == cantbl['name'])
				# score += np.asscalar(cantbl['score'][indx])
		elif ',' not in src:#	SINGLE SOURCES
			objlist.append(src)
			# indx = np.where(src == cantbl['name'])
			# score += np.asscalar(cantbl['score'][indx])
	set_objlist = list(set(objlist))
	objslist.append(objlist)
	# scorelist.append(score)
	score = 0.0
	set_score = 0.0
	for obj in objlist:
		indx = np.where(obj == cantbl['name'])
		score += np.asscalar(cantbl['score'][indx])
	for obj in set_objlist:
		indx = np.where(obj == cantbl['name'])
		set_score += np.asscalar(cantbl['score'][indx])
	f.write('{} {} {} {} {}\n'.format(obs, len(objlist), score, len(set_objlist), set_score))

allobjlist = []
for objs in objslist:
	for obj in objs:
		allobjlist.append(obj)
set_allobjlist = list(set(allobjlist))
allscore = countscore(allobjlist, cantbl)
set_allscore = countscore(set_allobjlist, cantbl)
f.write('{} {} {} {} {}\n'.format('TOTAL', len(allobjlist), allscore, len(set_allobjlist), set_allscore))
f.close()

restbl = ascii.read(path_save+'/result4initial.dat')
#------------------------------------------------------------
#	PHOTOMETRY + SCORE TABLE
#------------------------------------------------------------
scores = []
for src in photbl['sources']:
	score = 0.0
	if src == 'None':
		print('NO SOURCE', score)
		pass
	elif ',' in src:	#	MULTIPLE SOURCES
		for sr in src.split(','):
			indx = np.where(sr == cantbl['name'])
			score += np.asscalar(cantbl['score'][indx])
			print('MULTIPLE SOURCES', score)
	elif ',' not in src:#	SINGLE SOURCES
		indx = np.where(src == cantbl['name'])
		score += np.asscalar(cantbl['score'][indx])
		print('SINGLE SOURCE', score)
	scores.append(score)
col_scores = Column(np.array(scores), name='score')
photbl.add_column(col_scores)
photbl.write(path_save+'/phot+score4initial.dat', format='ascii', overwrite=True)