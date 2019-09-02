#	OBS. TABLE -> MEDIAN DEPTH+PHASE TABLE
#	2019.08.23	CREATED BY Gregory S.H. Paek
#============================================================
import numpy as np
from astropy.table import Table, vstack
from astropy.io import ascii                                                                                                          
import matplotlib.pyplot as plt
#============================================================
#	FUNCTION
#------------------------------------------------------------
def func(intbl, colist):
	newtbl = intbl[colist]
	newtbl = newtbl[np.argsort(newtbl['jd'])]
	return newtbl
#------------------------------------------------------------
def func2(intbl, trange, colist):
	intbl0 = func(intbl, colist)
	indx = np.where((intbl0['Phase']>=trange[0])&(intbl0['Phase']<=trange[1]))
	newtbl = intbl0[indx]
	phase_med = round(np.median(newtbl['Phase']), 3)
	ul_med = round(np.median(newtbl['ul']), 3)
	ul_std = round(np.std(newtbl['ul']), 3)
	outbl = Table(	[[intbl['obs'][0]], [intbl['filter'][0]], [phase_med], [ul_med], [ul_std]],
					names=('obs', 'filter', 'Phase', 'ul', 'ulstd'))
	return outbl
#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_save = '.'
path_base = '/data1/S190425z/1.result/table'
#------------------------------------------------------------
path_loao = path_base+'/obs_loao.dat'
path_lsgt = path_base+'/obs_lsgt.dat'
path_sao = path_base+'/obs_sao.dat'
path_squean = path_base+'/obs_squean.dat'
path_ukirt = path_base+'/obs_ukirt.dat'
path_kmtnet = path_base+'/obs_kmtnet.dat'
#------------------------------------------------------------
loaotbl = ascii.read(path_loao)
lsgtbl = ascii.read(path_lsgt)
saotbl = ascii.read(path_sao)
squeantbl = ascii.read(path_squean)
kmtnetbl = ascii.read(path_kmtnet)
ukirtbl = ascii.read(path_ukirt)
ukirtjtbl = ukirtbl[ukirtbl['filter']=='J']
ukirtktbl = ukirtbl[ukirtbl['filter']=='K']
#------------------------------------------------------------
colist = [ 'obs',
 'ra',
 'dec',
 'date-obs',
 'jd',
 'filter',
 'ul',
 'sources',
 'Phase']
#------------------------------------------------------------
# plt.hist(loaotbl['Phase'], bins='auto')
# plt.hist(lsgtbl['Phase'], bins='auto')
# plt.hist(saotbl['Phase'], bins='auto')
# plt.hist(squeantbl['Phase'], bins=np.arange(0, 5, 0.1)) 
# plt.hist(kmtnetbl['Phase'], bins=np.arange(0, 5, 0.1)) 
# plt.hist(ukirktbl['Phase'], bins=np.arange(0, 5, 0.1)) 
tblist = []
tblist.append(func2(loaotbl, [0, 0.6], colist))
tblist.append(func2(loaotbl, [0.9, 1.5], colist))
tblist.append(func2(loaotbl, [1.7, 2.5], colist))
tblist.append(func2(loaotbl, [2.5, 3.5], colist))
#------------------------------------------------------------
tblist.append(func2(lsgtbl, [0, 0.6], colist))
tblist.append(func2(lsgtbl, [1.0, 1.5], colist))
#------------------------------------------------------------
tblist.append(func2(saotbl, [1.0, 1.5], colist))
#------------------------------------------------------------
tblist.append(func2(squeantbl, [0, 0.5], colist))
tblist.append(func2(squeantbl, [0.5, 1.5], colist))
tblist.append(func2(squeantbl, [1.5, 2.5], colist))
tblist.append(func2(squeantbl, [2.5, 3.5], colist))
#------------------------------------------------------------
tblist.append(func2(kmtnetbl, [0, 1.0], colist))
tblist.append(func2(kmtnetbl, [1.0, 1.9], colist))
tblist.append(func2(kmtnetbl, [1.9, 3.0], colist))
#------------------------------------------------------------
tblist.append(func2(ukirtktbl, [0, 1.5], colist))
tblist.append(func2(ukirtjtbl, [1.5, 2.5], colist))
#------------------------------------------------------------
alltbl = vstack(tblist)
alltbl.write(path_save+'/obs_all.dat', format='ascii.fixed_width_two_line', overwrite=True)


'''
tftbl = ascii.read('/data1/S190425z/1.result/table/test/tableformat.txt')
formats = np.copy(tftbl['Format'])

for fmt in formats:
	try:
		alltbl.write('./test_table_{}.dat'.format(fmt), format=fmt, overwrite=True)
	except:
		print('fail {}'.format(fmt))
'''