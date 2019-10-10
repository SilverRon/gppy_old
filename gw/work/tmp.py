combimlist = []
for inim in imlist:
	combimlist.append(inim[:-7])
combimlist = list(set(combimlist))
combimlist.sort()
failist = []
for inim in combimlist:
	if inim not in photbl['obj']:
		failist.append(inim)
print(failist)
#----------
stdlist = []
for obj in combimlist: 
	stdlist.append(np.std(photbl[photbl['obj']==obj]['ul'])) 
plt.hist(stdlist, bins='auto')
#----------
obslist = []
for inim in photbl['image']:
	hdr = fits.getheader(inim)
	obslist.append(hdr['OBSERVAT'])
obslist = np.array(obslist)

newtbl = Table()
for key in photbl.keys():
	if key != 'obs':
		newtbl[key] = photbl[key]
	else:
		newtbl[key] = obslist


newtbl.write('phot-{}.dat'.format(input('name key:')), format='ascii')

'''
image obs obj ra dec date-obs jd filter stdnumb zp zper seeing skyval skysig ul mag magerr


'''


t0 = 2458598.84589
for intbl in [ssotbl, saaotbl, ctiotbl]:
	intbl['Phase'] = intbl['jd']-t0


#------------------------------------------------------------
intbl = ssotbl[ (ssotbl['Phase']<0.5) & (ssotbl['Phase']>0.0) ]
intbl = ssotbl[ (ssotbl['Phase']<1.5) & (ssotbl['Phase']>0.5) ]
intbl = ssotbl[ (ssotbl['Phase']<3.0) & (ssotbl['Phase']>1.5) ]
#------------------------------------------------------------
# intbl = saaotbl[ (saaotbl['Phase']<0.4) & (saaotbl['Phase']>0.0) ]
# intbl = saaotbl[ (saaotbl['Phase']<1.0) & (saaotbl['Phase']>0.4) ]
intbl = saaotbl[ (saaotbl['Phase']<1.0) & (saaotbl['Phase']>0.0) ]
intbl = saaotbl[ (saaotbl['Phase']<2.0) & (saaotbl['Phase']>1.0) ]
intbl = saaotbl[ (saaotbl['Phase']<3.0) & (saaotbl['Phase']>2.0) ]
#------------------------------------------------------------
intbl = ctiotbl[ (ctiotbl['Phase']<1.0) & (ctiotbl['Phase']>0.0) ]
intbl = ctiotbl[ (ctiotbl['Phase']<3.0) & (ctiotbl['Phase']>1.0) ]
#------------------------------------------------------------
intbl = loaotbl[ (loaotbl['Phase']<0.5) & (loaotbl['Phase']>0.0) ]
intbl = loaotbl[ (loaotbl['Phase']<1.5) & (loaotbl['Phase']>0.5) ]
intbl = loaotbl[ (loaotbl['Phase']<2.5) & (loaotbl['Phase']>1.5) ]
intbl = loaotbl[ (loaotbl['Phase']<3.5) & (loaotbl['Phase']>2.5) ]
#------------------------------------------------------------
intbl = lsgtbl[ (lsgtbl['Phase']<1.0) & (lsgtbl['Phase']>0.0) ]
intbl = lsgtbl[ (lsgtbl['Phase']<2.0) & (lsgtbl['Phase']>1.0) ]
#------------------------------------------------------------
intbl = sqtbl[ (sqtbl['Phase']<0.5) & (sqtbl['Phase']>0.0) ]
intbl = sqtbl[ (sqtbl['Phase']<1.5) & (sqtbl['Phase']>0.5) ]
intbl = sqtbl[ (sqtbl['Phase']<2.5) & (sqtbl['Phase']>1.5) ]
intbl = sqtbl[ (sqtbl['Phase']<3.5) & (sqtbl['Phase']>2.5) ]
#------------------------------------------------------------
intbl = utbl[ (utbl['Phase']<0.298) & (utbl['Phase']>0.0) ]
intbl = utbl[ (utbl['Phase']<1.3) & (utbl['Phase']>0.3) ]
intbl = utbl[ (utbl['Phase']<2.2) & (utbl['Phase']>1.9) ]
intbl = utbl[ (utbl['Phase']<6.4) & (utbl['Phase']>5.9) ]
#------------------------------------------------------------
def out(intbl):
	phasemed, ulmed = np.median( intbl['Phase'] ), np.median( intbl['ul'] )
	xer1 = phasemed - np.min( intbl['Phase'] )
	xer2 = np.max( intbl['Phase'] ) - phasemed
	yer1 = ulmed - np.min( intbl['ul'] )
	yer2 = np.max( intbl['ul'] ) - ulmed

	phasemed = round(phasemed, 3)
	ulmed = round(ulmed, 3)
	xer1 = round(xer1, 3)
	xer2 = round(xer2, 3)
	yer1 = round(yer1, 3)
	yer2 = round(yer2, 3)

	print(phasemed, ulmed, xer1, xer2, yer1, yer2)

out(intbl)