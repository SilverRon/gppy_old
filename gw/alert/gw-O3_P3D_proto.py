#   UTILIZE GW P_3D LOCALIZATION MAP
#   BASED ON https://github.com/lpsinger/gw-galaxies/blob/master/gw-galaxies.ipynb
#   2019.07.19 MADE		BY Gregory S.H. Paek
#	2019.XX.XX UPDATED	BY Gregory S.H. Paek
#============================================================
#	MODULE
#------------------------------------------------------------
from ligo.skymap.postprocess import find_greedy_credible_levels
import gcn
import gcn.handlers
import gcn.notice_types
import healpy as hp
import numpy as np
#------------------------------------------------------------
import astropy.utils.data
from matplotlib import pyplot as plt
import numpy as np
import healpy as hp
#------------------------------------------------------------
#	FOR THE GALAXY CROSS MATCHING
from astropy.table import Table, vstack, hstack, Column
import astropy.units as u
from astropy.coordinates import SkyCoord
import ligo.skymap.plot
from scipy.stats import norm
import scipy.stats
#============================================================
url			= 'https://dcc.ligo.org/public/0146/G1701985/001/bayestar.fits.gz'
# filename	= astropy.utils.data.download_file(url)
filename	= 'S190425z-Update-bayestar.fits.gz'
hpx, hdr	= hp.read_map(filename, verbose=True, h=True)
hdr			= dict(hdr)
prob, distmu, distsigma, distnorm = hp.read_map(filename,
												#field=[0, 1, 2, 3],
												field=range(4),
												dtype=('f8', 'f8', 'f8', 'f8'))
credible_levels = find_greedy_credible_levels(prob)
indx = np.where(credible_levels <= 0.9)
# prob, distmu,  distsigma, distnorm = prob[indx], distmu[indx],  distsigma[indx], distnorm[indx]


npix		= len(prob)
nside		= hp.npix2nside(npix)
pixarea		= hp.nside2pixarea(nside)		
pixarea_deg2= hp.nside2pixarea(nside, degrees=True)
#------------------------------------------------------------
path_cat	= '/home/sonic/Research/cat/GLADE2.3/GLADE_2.3+2MASS_PSC+identi_name.dat'
gldtbl0		= Table.read(path_cat, format='ascii')
n=1
gldtbl		= gldtbl0[	(gldtbl0['dist']<=hdr['DISTMEAN']+n*hdr['DISTSTD'])&
						(gldtbl0['dist']>=hdr['DISTMEAN']-n*hdr['DISTSTD'])]
gldcoord	= SkyCoord(ra=gldtbl['ra']*u.deg, dec=gldtbl['dec']*u.deg)
ngld		= np.size(gldtbl)

probdencol	= Column(np.zeros(ngld, dtype='f4'), name='dP_dV')
probcol		= Column(np.zeros(ngld, dtype='f4'), name='P')
probdenAcol	= Column(np.zeros(ngld, dtype='f4'), name='dP_dA')
probAcol	= Column(np.zeros(ngld, dtype='f4'), name='P_A')
gldtbl.add_columns([probdencol, probcol, probdenAcol, probAcol])
#------------------------------------------------------------
#get coord of max prob density for plotting purposes
ipix_max	= np.argmax(prob)
theta_max, phi_max = hp.pix2ang(nside, ipix_max)
ra_max, dec_max	= np.rad2deg(phi_max), np.rad2deg(0.5 * np.pi - theta_max)
# ra_max, dec_max	= hp.pix2ang(nside, ipix_max, lonlat=True)
center		= SkyCoord(ra=ra_max*u.deg, dec=dec_max*u.deg)
# ipix		= np.arange(0, npix, step=1)
# theta, phi	= hp.pix2ang(nside, ipix)
# ra, dec		= np.rad2deg(phi), np.rad2deg(0.5 * np.pi - theta)

'''
#------------------------------------------------------------
#	PLOT
#------------------------------------------------------------
#plot up the sky localization and overplot the galaxies
#give the students most of this. Later they will have to add to it in order to overplot the galaxies
ax = plt.axes(	[0.05, 0.05, 0.9, 0.9],
				projection='astro globe',
				center=center)
ax.grid()
ax.imshow_hpx(filename, cmap='cylon')
ax.plot(center.ra.deg, center.dec.deg,
		transform=ax.get_transform('world'),
		marker=ligo.skymap.plot.reticle(inner=0, outer=1),
		markersize=10, markeredgewidth=2)
'''
#------------------------------------------------------------
#calc hp index for each galaxy
theta	= 0.5 * np.pi - gldcoord.dec.to('rad').value
phi		= gldcoord.ra.to('rad').value
ipix	= hp.ang2pix(nside, theta, phi)
cumP2D	= np.cumsum(prob[np.argsort(-1*prob)])[ipix]
'''
ipix[cumP2D < 0.9]
ipix[cumP2D < 0.9]
dp_dV	= prob[ipix[cumP2D < 0.9]] * distnorm[ipix[cumP2D < 0.9]] * norm(distmu[ipix[cumP2D < 0.9]],distsigma[ipix[cumP2D < 0.9]]).pdf(gldtbl[cumP2D < 0.9]['dist'])/pixarea 
'''
#calc probability (P_2D)
dp_dA		= prob[ipix]/pixarea
dp_dA_deg2	= prob[ipix]/pixarea_deg2
#calc probability density per volume for each galaxy
dp_dV	= prob[ipix] * distnorm[ipix] * norm(distmu[ipix],distsigma[ipix]).pdf(gldtbl['dist'])/pixarea 
#------------------------------------------------------------
gldtbl['dP_dV']	= dp_dV
gldtbl['dP_dA']	= dp_dA
gldtbl['Prob']	= gldtbl['dist']**2 * 10**(-0.4*gldtbl['K']) * gldtbl['dP_dV']
#------------------------------------------------------------
cantbl			= gldtbl[	(gldtbl['K']!=-99.0)&
							(gldtbl['dist']!=-99.0)&
							(gldtbl['Prob']!=0.0)]
cantbl['Prob']	= cantbl['Prob']/np.sum(cantbl['Prob'])
cantbl.sort('Prob')
cantbl.reverse()

subtbl			= cantbl[np.cumsum(cantbl['Prob'])<=0.9]
# cantbl['']		= cantbl[np.argsort(-1*cantbl['Prob'])]
