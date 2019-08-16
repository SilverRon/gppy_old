#	KILONOVA MODEL SPECTRUM TO K-CORRECTED MAGNITUDE BASED ON Kasen+2017
#	REFERENCE	: http://www.astro.ljmu.ac.uk/~ikb/research/mags-fluxes/
#------------------------------------------------------------
#	2019.08.??	BASED ON	Joonho Kim
#	2019.08.14	MODIFIED BY Gregory S.H.Paek
#============================================================
import h5py
import numpy as np
import bisect
import matplotlib.pyplot as plt 
import os,sys
import glob
from scipy.integrate import simps
from astropy.io import ascii
from astropy import units as u
from astropy import constants as const
from speclite import filters
from astropy.table import Table, Row, Column, vstack
#============================================================
#	FUNCTION
#------------------------------------------------------------
def model2mag(mdl, filters, z=0, h=0.7, path_save='.'):
	print (mdl)
	#------------------------------------------------------------
	H0 = 100*h								#	[km/s/Mpc]
	d = 1e6*const.c.to_value('km/s')*z/H0	#	[Mpc]*1e6 -> [pc] 
	#------------------------------------------------------------
	tblist = []
	outname = path_save+'/{0}_z{1}_mag.dat'.format(os.path.basename(mdl), z)
	fin = h5py.File(mdl, 'r')
	'''
	fin.keys()
	<KeysViewHDF5 ['Lnu', 'mu', 'nu', 'time']>
	'''
	nu = np.array(fin['nu'], dtype='d')			#	[Hz]
	times = np.array(fin['time'])/3600.0/24.0	#	[sec] -> [day]
	Lnu_all = np.array(fin['Lnu'], dtype='d')	#	SPECIFIC LUMINOSITY	[erg/s/Hz]

	for t in times:
		onetbl = Table()
		# onetbl.add_column([t], names=('t'))
		onetbl.add_column(Column(name='t', data=[t]))
		# magarr = np.zeros(len(filters))
		it = bisect.bisect(times, t)
		Lnu = Lnu_all[it-1,:]
		#	CHANGE VARAIABLE : NU -> LAMDA
		lam = (1+z)*const.c.to_value('Angstrom/s')/nu		#	[Ang]
		Llam = Lnu*nu**2.0/const.c.to_value('Angstrom/s')	#	SPECIFIC LUMINOSITY [erg/s/Ang]
		indx_sort = np.argsort(lam)
		lam = lam[indx_sort]
		Llam = Llam[indx_sort]
		Lflux = Llam/(4*np.pi*(const.pc.to_value('cm'))**2)
		#	F_MU TO AB MAGNITUDE FOR EACH FILTERS
		for j, filte in enumerate(filters):
			filtbl = ascii.read(filte)
			lamda = filtbl['lamda']
			trans = filtbl['transmission']
			
			filte_intp = np.interp(lam, lamda, trans)
			intensity1 = simps(Lflux*filte_intp*lam, lam)
			intensity2 = simps(filte_intp/lam, lam)
			fmu = intensity1/intensity2/const.c.to_value('Angstrom/s')
			mag = -2.5*np.log10(fmu)-48.6
			
			if (mag == float('inf')) | (mag == float('nan')):
				appmag = -99.0
			else:
				appmag = abs2app(mag, d)
			# if (mag == float('inf')) | (mag == float('nan')): mag = -99
			onetbl.add_column(Column(name=os.path.basename(filte)[0], data=[appmag]))
		tblist.append(onetbl)
	lctbl = vstack(tblist)
	if path_save == False:
		pass
	else:
		lctbl.write(outname, format='ascii', overwrite=True)
	return lctbl
#------------------------------------------------------------
def abs2app(absmag, d):
	'''
	absmag	:	ABSOLUTE MAGNITUDE	[mag]
	d		:	DISTANCE			[pc]
	'''
	appmag = absmag+5*np.log10(d/10)
	return appmag
#============================================================
#	PATH
#------------------------------------------------------------
path_models = '/home/sonic/Research/yourpy/gppy/table/Kasen_Kilonova_Models_2017/kilonova_models'
path_filters = '/home/sonic/Research/yourpy/gppy/table/filter_transmission'
path_save = '/data1/S190425z/1.result/table'
#------------------------------------------------------------
#	GW170817	(http://simbad.u-strasbg.fr/simbad/sim-id?Ident=NGC+4993)
z = 0.009787
models = glob.glob(path_models+'/*h5')
filters = glob.glob(path_filters+'/*.dat')

# mdl = models[0]
plt.close('all')
for i, mdl in enumerate(models):
	print('[{}/{}]'.format(i+1, len(models)))
	# if ('m0.025' in mdl) & ('vk0.3' in mdl) & ('Xlan1e-4.0' in mdl):
	if ('m0.04' in mdl) & ('vk0.15' in mdl) & ('Xlan1e-1.5' in mdl):
		outname = path_save+'/{0}_z{1}_mag.dat'.format(os.path.basename(mdl), z)
		# if outname not in glob.glob(path_save+'/*'):
		lctbl = model2mag(mdl=mdl, filters=filters,
						z=z, h=0.7, path_save=path_save)
		# else:
			# lctbl = ascii.read(outname)
		rlctbl = lctbl[lctbl['i'] != -99]
		# plt.scatter(rlctbl['t'], rlctbl['r'])
		plt.plot(rlctbl['t'], rlctbl['i'], label='{}'.format(os.path.basename(mdl)))
# plt.ylim([22, 15])



droutbl = ascii.read(path_save+'/lc_gw170817_Drout.dat')
rdroutbl = droutbl[droutbl['filter']=='i']
plt.scatter(rdroutbl['delmjd'], rdroutbl['mag'], marker='o', c='grey', label='Drout+17')
plt.xlim([-1, 11])
plt.ylim([25, 10])

plt.minorticks_on()
plt.tight_layout()
