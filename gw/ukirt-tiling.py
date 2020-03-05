#   UKIRT TILEING
#   2019.08.17	CREATED BY Gregory S.H. Paek
#============================================================    
#   MODULE
import os, sys, glob
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib
from astropy.table import Table, Column, MaskedColumn, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier 
from astropy.coordinates import Angle
import astropy.coordinates as coord
#============================================================
#	FUNCTION
#------------------------------------------------------------
def maketile(ra, dec, fov):
	'''
				   1*----------*4
				   *            *
				 2*--------------*3
				
				RA -->
	'''
	dec_upper = dec + fov/2.
	dec_lower = dec - fov/2.
	ra1 = ra - (fov/2.)/np.cos((dec_upper*u.deg).to_value('radian'))
	ra2 = ra - (fov/2.)/np.cos((dec_lower*u.deg).to_value('radian'))
	ra3 = ra + (fov/2.)/np.cos((dec_lower*u.deg).to_value('radian'))
	ra4 = ra + (fov/2.)/np.cos((dec_upper*u.deg).to_value('radian'))
	tile = dict(center=dict(ra=ra,dec=dec),
				pt1=dict(ra=ra1,dec=dec_upper),
				pt2=dict(ra=ra2,dec=dec_lower),
				pt3=dict(ra=ra3,dec=dec_lower),
				pt4=dict(ra=ra4,dec=dec_upper),
				score=0.0)
	return tile
#------------------------------------------------------------
def totalscore(tile, intbl, scorekey):
	'''
	APPROXIMATE SCORE (SIMPLIFYING TILTED EDGE)
	'''
	ramean12 = np.mean([tile['pt1']['ra'], tile['pt2']['ra']])
	ramean34 = np.mean([tile['pt3']['ra'], tile['pt4']['ra']])
	dec_upper = tile['pt1']['dec']
	dec_lower = tile['pt2']['dec']
	indx = np.where((intbl['ra']>ramean12)&
					(intbl['ra']<ramean34)&
					(intbl['dec']>dec_lower)&
					(intbl['dec']<dec_upper))
	totscore = np.sum(intbl[scorekey][indx])
	tile['score'] = totscore
	return tile
#------------------------------------------------------------
def ratiling(tilecent, fov, intbl, scorekey='P_2D', iter=100, cut=0.0):
	racent, decent = tilecent['center']['ra'], tilecent['center']['dec']
	step = fov/np.cos((decent*u.deg).to_value('radian'))
	tilelist = []
	for i in np.arange(-1*iter, iter+1, 1):
		tile = maketile(racent+i*step, decent, fov)
		tile = totalscore(tile, intbl, scorekey)
		if tile['score'] > cut:
			tilelist.append(tile)
	return tilelist
#------------------------------------------------------------
def drawtile(tile, ax):
	ra1, ra2, ra3, ra4 = tile['pt1']['ra'], tile['pt2']['ra'], tile['pt3']['ra'], tile['pt4']['ra']
	dec1, dec2, dec3, dec4 = tile['pt1']['dec'], tile['pt2']['dec'], tile['pt3']['dec'], tile['pt4']['dec']
	ax.plot([ra1, ra2, ra3, ra4, ra1], [dec1, dec2, dec3, dec4, dec1], c='dodgerblue')
#------------------------------------------------------------
def drawtile4table(racent, decent, fov, ax):
	tile = maketile(racent, decent, fov)
	drawtile(tile, ax)
#------------------------------------------------------------
def pickopt(raoslist):
	scores = []
	for i in range(len(raoslist)):
		totscore = 0.0
		tilenumb = 0
		offset = raoslist[i]['offset']
		tilenumb += len(raoslist[i]['tiles'])
		for t in raoslist[i]['tiles']:
			totscore += t['score']
		try:
			scores.append(1/((totscore)**2)/tilenumb)
		except:
			scores.append(0.0)
	indx_opt = np.where(scores == np.max(scores))
	opt_tiles = raoslist[indx_opt[0][0]]
	return opt_tiles
#------------------------------------------------------------
def tiletable(tilelists):
	tblist = []
	ralist = []
	delist = []
	sclist = []
	ct = 0
	for i, tilelist in enumerate(tilelists):
		for j, tile in enumerate(tilelist['tiles']):
			ralist.append(round(tile['center']['ra'], 3))
			delist.append(round(tile['center']['dec'], 3))
			sclist.append(tile['score'])
			ct += 1
	tiletbl = Table([np.arange(1, ct+1, 1), ralist, delist, sclist],
					names=('numb', 'ra', 'dec', 'score'))
	tiletbl = tiletbl[np.argsort(-1*tiletbl['score'])]
	tiletbl['numb'] = np.arange(1, len(tiletbl)+1, 1)
	return tiletbl
#------------------------------------------------------------
def coord2msb(c):
	cstr = c.to_string('hmsdms')
	cstr = cstr.replace('h', ' ')
	cstr = cstr.replace('m', ' ')
	cstr = cstr.replace('s', ' ')
	cstr = cstr.replace('d', ' ')
	cstr = cstr.replace(' -', '- ')
	cstr = cstr.replace(' +', '+ ')
	return cstr
#============================================================
#	UKIRT ONE TILE -> 4 TILES FOR MSB
def fourtiles(ra, dec, fov):
	'''
	1       4
	*-------*
	|       |
	|   x   |
	|       |
	*-------*
	2       3

	RA ->

	'''
	dec_upper = dec+fov/8.
	dec_lower = dec-fov/8.
	ra1 = ra - (fov/8.)/np.cos((dec_upper*u.deg).to_value('radian'))
	ra2 = ra - (fov/8.)/np.cos((dec_lower*u.deg).to_value('radian'))
	ra3 = ra + (fov/8.)/np.cos((dec_lower*u.deg).to_value('radian'))
	ra4 = ra + (fov/8.)/np.cos((dec_upper*u.deg).to_value('radian'))

	ralist = [ra1, ra2, ra3, ra4]
	delist = [dec_upper, dec_lower, dec_lower, dec_upper]

	return ralist, delist
#------------------------------------------------------------
def fourtiles4ukirt(ra, dec):

	a = 13.65/60.
	b = 12.83/60.
	c = 4.25/60.
	r = c/2
	step = (a+b)/2

	dec_up = dec+step/2.
	dec_lo = dec-step/2.
	ra1 = ra - (step/2.)/np.cos((dec_up*u.deg).to_value('radian'))
	ra2 = ra - (step/2.)/np.cos((dec_lo*u.deg).to_value('radian'))
	ra3 = ra + (step/2.)/np.cos((dec_lo*u.deg).to_value('radian'))
	ra4 = ra + (step/2.)/np.cos((dec_up*u.deg).to_value('radian'))

	ralist = [ra1, ra2, ra3, ra4]
	delist = [dec_up, dec_lo, dec_lo, dec_up]

	return ralist, delist
#------------------------------------------------------------
def tbl2msb(comtbl, eventname, path_save):
	f = open(path_save+'/'+eventname+'-UKIRT-tiling4survey.dat', 'a')
	for l in range(len(comtbl)):
		# name = np.asscalar(comtbl['tile'][l])
		name = comtbl['tile'][l].item()
		c = SkyCoord(ra=comtbl['ra'][l]*u.degree, dec=comtbl['dec'][l]*u.degree, frame='icrs')
		cstr = c.to_string('hmsdms')
		cstr = cstr.replace('h', ':')
		cstr = cstr.replace('m', ':')
		cstr = cstr.replace('s', '')
		cstr = cstr.replace('d', ':')
		cstr = cstr.replace('+', '+')
		cstr = cstr.replace('-', '-')
		line = 'BASE {} {} RJ\n'.format(name, cstr)
		f.write(line)
		print(line)
		#	GUIDE STAR CHECK
		ra, dec = comtbl['ra'][l], comtbl['dec'][l]
		respond = guidering4ukirt(ra, dec)
		try:
			gtbl = respond[:2]
			for n in range(len(gtbl)):
				part = name.split('-'); part[0] = 'GD{}'.format(n+1)
				gname = '-'.join(part)
				gra, gdec = gtbl['RAJ2000'][n], gtbl['DEJ2000'][n]
				gc = SkyCoord(gra, gdec, frame='icrs', unit='deg')
				cstr = gc.to_string('hmsdms')
				cstr = cstr.replace('h', ':')
				cstr = cstr.replace('m', ':')
				cstr = cstr.replace('s', '')
				cstr = cstr.replace('d', ':')
				cstr = cstr.replace('+', '+')
				cstr = cstr.replace('-', '-')
				if (n+1) == 1:
					gnumb = ''
				elif (n+1) == 2:
					gnumb = '2'
				line = 'GUIDE{} {} {} {} FK5\n'.format(gnumb, gname, gra, gdec)
				print(line)
				f.write(line)
		except:
			pass

	f.close()
#------------------------------------------------------------
def guidering4ukirt(ra, dec):
	r = 4.25/2.
	# r = 10
	Vizier.ROW_LIMIT    = -1
	query   = Vizier.query_region(SkyCoord(ra=ra, dec=dec,
								unit=(u.deg, u.deg), frame='icrs'),
								width=str(r)+'m', catalog=["APASS9"])
	if len(query) != 0:
		intbl = query[0]
		gtbl = intbl[intbl['r_mag']<13.0]
		if len(gtbl) != 0:
			guidetbl = gtbl[np.argsort(gtbl['r_mag'])]
			return guidetbl
		else:
			return False
	else:
		return False

#------------------------------------------------------------
'''

<--------->
a   b  
----   ----
|  |   |  |
----   ----
     X ---> c 
----   ----
|  |   |  |
----   ----

==========> RA

'''
a = 13.65
b = 12.83
c = 4.25
r = c/2
step = (a+b)/2

fov = (5*a+3*b)/120.		#	[arcmin]

#============================================================
#	SETTING
#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
# os.system('locate */S*-*healpix*.dat')
os.system('ls *.dat')
path_healpix = input('HEALPIX TABLE\t: ')
# path_healpix = 'S190814bv-Update_healpix_70.dat'
eventname = input('EVENTNAME (S190814bv)\t: ')
if eventname == '': eventname = 'S190814bv'
path_save = input('SAVE PATH (.)\t: ')
if path_save == '': path_save = '.'
#------------------------------------------------------------
# fov = 0.85	#	UKIRT TILE(4 PAWS) ~ 0.85 [deg]
htbl = ascii.read(path_healpix)
# htbl = htbl[(htbl['dec']>0)&(htbl['dec']<60)]
indx_max = np.where(htbl['P_2D']==np.max(htbl['P_2D']))
# ramax, demax = np.asscalar(htbl['ra'][indx_max]), np.asscalar(htbl['dec'][indx_max])
ramax, demax = htbl['ra'][indx_max].item(), htbl['dec'][indx_max].item()
#------------------------------------------------------------
#	Declination ITERATION
#------------------------------------------------------------
tilelists = []
# for decent in np.arange(np.min(htbl['dec'])-fov*0.5, np.max(htbl['dec'])+fov*1.5, fov):
for i, decent in enumerate(np.arange(np.min(htbl['dec'])-fov*2.0, np.max(htbl['dec'])+fov*2.0, fov)):
	print('PROCESS [{}/{}]'.format(i+1, len(np.arange(np.min(htbl['dec'])-fov*0.5, np.max(htbl['dec'])+fov*1.5, fov))))
	raoslist = []	#	RA OFFSET LIST
	tilecent = maketile(ramax, decent, fov)
	#------------------------------------------------------------
	#	Right Accession ITERATION (FINDING OPTIMIZED OFFSET)
	#------------------------------------------------------------
	for offset in np.arange(-0.5, 0.5+fov/100, fov/100):
		newtilecent = tilecent
		newtilecent['center']['ra'] += offset
		newtilecent = totalscore(newtilecent, htbl, 'P_2D')
		ratilelist = ratiling(newtilecent, fov, htbl, cut=0.0)
		raoslist.append(dict(offset=offset, tiles=ratilelist))
		opt_tiles = pickopt(raoslist)
	if (len(opt_tiles['tiles']) != 0) & (newtilecent['score'] != 0.0):
		tilelists.append(opt_tiles)
#------------------------------------------------------------
#	TILE TABLE & MSB FILE
#------------------------------------------------------------
tiletbl = tiletable(tilelists)
onesig = np.median(tiletbl['score']) - np.std(tiletbl['score'])
tiletbl.write(path_save+'/{}-UKIRT-tiling.dat'.format(eventname), format='ascii', overwrite=True)
'''
os.system('rm {}/*4script.txt'.format(path_save))
f = open(path_save+'/'+eventname+'-UKIRT-tiling-4script.txt', 'a')
for i in range(len(tiletbl)):
	if len(str(i+1)) < 3:
		numbering = str(i+1)
		for j in range(3-len(str(i+1))):
			numbering = '0'+numbering
	else:
		numbering = i+1
	tname = '{}-{}'.format(eventname, numbering)
	c = SkyCoord(tiletbl['ra'][i], tiletbl['dec'][i], frame='icrs', unit='deg')
	cstr = coord2msb(c)
	line = '{} {} RJ\n'.format(tname, cstr)
	print(line)
	f.write(line)
f.close()
'''
#------------------------------------------------------------
#	MSB FILE FOR SURVEY CONTAINER
#------------------------------------------------------------
tblist = []
for i in range(len(tiletbl)):
	n = i+1
	newname = eventname+'-'
	if len(str(n)) < 3:
		for j in range(3-len(str(n))):
			newname = newname+'0'
	newname = newname+str(n)+'-'
	ra, dec = tiletbl['ra'][i], tiletbl['dec'][i]
	ralist, delist = fourtiles4ukirt(ra, dec)

	namelist = []
	for k in ['1', '2', '3', '4']:
		namelist.append(newname+k)
	septbl = Table([namelist, ralist, delist],
					names=('tile', 'ra', 'dec'))
	tblist.append(septbl)
comtbl = vstack(tblist)
tbl2msb(comtbl, eventname, path_save)
#------------------------------------------------------------
#	GUIDE STAR SELTECTION WITH 2MASS CATALOG (PREPARE)
#------------------------------------------------------------
#============================================================
#	PLOT TILES ON HEALPIX SKYMAP
#============================================================
plt.close('all')
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
#------------------------------------------------------------
#	HEALpix
#------------------------------------------------------------
ax.scatter(htbl['ra'], htbl['dec'], marker='o', s=1, c=-1*htbl['P_2D'], cmap='OrRd')
#------------------------------------------------------------
#	TILE
#------------------------------------------------------------
for i in range(len(tiletbl)):
	numb = tiletbl['numb'][i]
	racent, decent = tiletbl['ra'][i], tiletbl['dec'][i]
	drawtile4table(racent, decent, fov, ax)
	ax.text(racent+0.05, decent-0.05, numb, fontsize=15)
#------------------------------------------------------------
#	PAW
#------------------------------------------------------------
ax.scatter(comtbl['ra'], comtbl['dec'], marker='x', c='k')
'''
for i, tiles in enumerate(tilelists):
	print('tilelists [{}/{}]'.format(i+1, len(tilelists)))
	for j, t in enumerate(tiles['tiles']):
		# print('tiles [{}/{}]'.format(j+1, len(tiles)))
		drawtile(t, ax)
	# nextcom = input('NEXT?:')
'''
plt.xlabel('RA [deg]', fontsize=20)
plt.ylabel('Dec [deg]', fontsize=20)
plt.minorticks_on()
plt.tight_layout()
plt.gca().invert_xaxis()
plt.axis(option='equal')
plt.savefig(path_save+'/{}-UKIRT-tiling-new.png'.format(eventname), overwrite=True)
