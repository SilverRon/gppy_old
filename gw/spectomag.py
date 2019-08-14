#Kasen+2017 kilonova model spectrum to K-corrected broadband magnitude
#ex python mkspec.py redshift

import h5py
import numpy as np
import bisect
import matplotlib.pyplot as plt 
import os,sys
import glob
from scipy.integrate import simps
from speclite import filters

z=0
z=sys.argv[1]
path = "./Kasen_Kilonova_Models_2017-master/kilonova_models/*.h5"
filelist = glob.glob(path)
c = 2.99792458e10
ca = 2.99792458e18
filters = ['u','g','r','i','z','U','B','V','R','I']

for nfilelist in filelist:
	print (nfilelist)
#	f = open(nfilelist+'_mag'+'.dat','w')
	f = open(nfilelist+'_z'+str(z)+'_mag'+'.dat','w')
	fin = h5py.File(nfilelist,'r')
	nu = np.array(fin['nu'],dtype='d')
	times = np.array(fin['time'])
	times = times/3600.0/24.0
	Lnu_all = np.array(fin['Lnu'],dtype='d')

	for time in times:
		magarr = np.zeros((10))
		it = bisect.bisect(times,time)
		Lnu = Lnu_all[it-1,:]
		lam  = (1+float(z))*(2.99e10)/nu*1e8
		Llam = Lnu*nu**2.0/(2.99e10)/1e8
		sort=np.argsort(lam)
		lam=lam[sort]
		Llam=Llam[sort]
		Lflux = Llam/(4*3.141592*(3.086e+19)**2)

		for j,nfilter in enumerate(filters):
			filterdatname = glob.glob('*'+nfilter+'.dat')
			filterdat = np.loadtxt(filterdatname[0])
			lamf = filterdat[:,0]
			filt = filterdat[:,1]
			filterint = np.interp(lam,lamf,filt)
			intensity1 = simps(Lflux*filterint*lam,lam)
			intensity2 = simps(filterint/lam,lam)
			flam = intensity1/intensity2/ca
			magarr[j] = -2.5*np.log10(flam)-48.6

		f.write('{:7.3f}	{:7.3f}	{:7.3f}	{:7.3f}	{:7.3f}	{:7.3f}	{:7.3f}	{:7.3f}	{:7.3f}	{:7.3f}	{:7.3f}'.format(time, magarr[0], magarr[1], magarr[2], magarr[3], magarr[4], magarr[5], magarr[6], magarr[7], magarr[8], magarr[9])+'\n')
	f.close

#plt.plot(lam,Llam)
#plt.ion()
#plt.show()

#		f = open(i+'_'+str(round(time,2))+".dat",'w')
#		f.write(str(round(time,2))+' '+str(round(time,2)))
#		for j in range(len(lam)):
#			f.write('{:e} {:e}'.format(lam[j], Llam[j])+'\n')
#		f.close


