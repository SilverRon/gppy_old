from astropy import units as u
from astropy.coordinates import SkyCoord


def deg2hmsdms(radeg, dedeg):
	c = SkyCoord(ra=radeg*u.degree, dec=dedeg*u.degree, frame='icrs')
	hmsdms = c.to_string('hmsdms')
	for i in ['h', 'm', 'd']:
		hmsdms = hmsdms.replace(i, ':')
	hmsdms = hmsdms.replace('s', '')
	print(hmsdms)

radeg, dedeg = 323.19, 4.53		#	[deg]
step = 1.05						#	[deg]

deg2hmsdms(radeg+step, dedeg+step)
deg2hmsdms(radeg+step, dedeg-step)
deg2hmsdms(radeg-step, dedeg+step)
deg2hmsdms(radeg-step, dedeg-step)
