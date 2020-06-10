#	rts_vis MAKER WITH astroplan MODULE
#	2019.12.18	CREATED BY	Gregory S.H. Paek
#============================================================
#	MODULES
#------------------------------------------------------------
import os, glob
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import warnings
warnings.filterwarnings(action='ignore')
#------------------------------------------------------------
from astropy.time import Time
from astropy.io import ascii
from astropy.coordinates import SkyCoord, get_moon
# from pytz import timezone
import pytz
#------------------------------------------------------------
from astroplan import Observer
from astropy.coordinates import EarthLocation
from astroplan.plots import plot_airmass
#============================================================
#	FUNCTION
#------------------------------------------------------------
def ut2local(observer, time):
	isoform = observer.astropy_time_to_datetime(time).isoformat()
	hhmm = (isoform.split('T')[1])[0:5]
	return hhmm
#------------------------------------------------------------
def lct2utc(lct):
	utctuple = lct.utctimetuple()
	utc = '{}-{}-{} {}:{}:{}'.format(utctuple[0], utctuple[1], utctuple[2], utctuple[3], utctuple[4], utctuple[5])
	return utc
#------------------------------------------------------------
def callobserver(obs, obstbl, description=''):
	indx_obs = np.where( obs == obstbl['name'] )
	longitude, latitude, elevation = obstbl[indx_obs]['longitude(E+)'], obstbl[indx_obs]['latitude(N+)'], obstbl[indx_obs]['altitude']
	tz = (obstbl[indx_obs]['timezone']).item()
	location = EarthLocation.from_geodetic(longitude, latitude, elevation)

	observer = Observer(name=obs,
						location=location,
						# pressure=0.615 * u.bar,
						# relative_humidity=0.11,
						# temperature=0 * u.deg_C,
						timezone=tz,
						description=description)
						# description="SAO 1-m Telescope on Seoul National University, Korea")
	return observer
#------------------------------------------------------------
def rts_vis_maker(c, observer, y, m, d, path_out='.'):
	#------------------------------------------------------------
	#	TIME
	#------------------------------------------------------------
	# time = Time('2020-01-01 12:00:00')		# UTC
	pst = pytz.timezone(observer.timezone.zone)
	lct_input = pst.localize(datetime(y, m, d, 1, 0, 0))	#	LOCAL TIME
	try:
		lct = pst.localize(datetime(y, m, d+1, 1, 0, 0))	#	LOCAL TIME
	except:
		lct = pst.localize(datetime(y, m+1, 1, 1, 0, 0))
	time = Time(lct2utc(lct))
	dts = pst.localize(time.datetime).dst().total_seconds()/3600.
	if dts >= 0:
		dts = '+{}'.format(dts)
	else:
		dts = str(dts)
	#------------------------------------------------------------
	#	WRITE rts_vis
	#------------------------------------------------------------
	f = open('{}/rts_vis_{}_IMSNG_{}.txt'.format(path_out, ''.join(lct_input.isoformat().split('T')[0].split('-')), obs), 'w')
	header = False
	#------------------------------------------------------------
	#	SUN INFO.
	#------------------------------------------------------------
	sunset_tonight = observer.sun_set_time(time, horizon=-18*u.deg, which='nearest')
	sunrise_tonight = observer.sun_rise_time(time, horizon=-18*u.deg, which='nearest')
	# night_time = sunset_tonight + ((sunrise_tonight-0.25) - (sunset_tonight+0.25))*np.linspace(0, 1, 1000)
	for i, target in enumerate(c):
		#------------------------------------------------------------
		#	ONE TARGET
		#------------------------------------------------------------
		try:
			c_rise = observer.target_rise_time(time, target, horizon=constraint_alt, which='previous')
			c_set = observer.target_set_time(time, target, horizon=constraint_alt, which='next')
			all_up_start = np.max([c_rise])
			all_up_end = np.min([c_set])
			start = np.max([sunset_tonight, all_up_start])
			end = np.min([sunrise_tonight, all_up_end])
			visible_time = start + (end - start)*np.linspace(0, 1, 100)
			#	POSITION
			altaz = observer.altaz(visible_time, target)
			alt = altaz.alt
			az = altaz.az
			night_time = sunset_tonight + ((sunrise_tonight) - (sunset_tonight))*np.linspace(0, 1, 1000)
			target_time = c_rise + (c_set - c_rise)*np.linspace(0, 1, 1000)
			altaz_target = observer.altaz(target_time, target)
			alt_target = altaz_target.alt
			az_target = altaz_target.az
			#------------------------------------------------------------
			#	MOON INFO.
			#------------------------------------------------------------
			moon_altaz = observer.moon_altaz(visible_time)
			moon_alt, moon_az = moon_altaz.alt, moon_altaz.az

			c_moon = SkyCoord(moon_az, moon_alt, frame='altaz')
			moon_illumi = round(np.max(observer.moon_illumination(night_time)), 2)
			c_target = SkyCoord(az, alt, frame='altaz')
			sep_moon = c_target.separation(c_moon)
			#------------------------------------------------------------
			# c_moon_radec = get_moon(visible_time, location=observer.location).transform_to('icrs')
			# sep_moon = target.separation(c_moon_radec)
			#------------------------------------------------------------
			#	ALT, MOON DIST. CONSTRAINT
			#------------------------------------------------------------
			if (np.max(sep_moon) >= constraint_moon) & (np.max(alt) >= constraint_alt):
				transit_ut = Time((altaz_target[alt_target == np.max(alt_target)]).obstime.iso.item())
				#------------------------------------------------------------
				name = intbl['obj'][i]
				ra, dec = intbl['ra'][i], intbl['dec'][i]
				#------------------------------------------------------------
				rise_local = ut2local(observer, c_rise)
				transit_local = ut2local(observer, transit_ut)
				set_local = ut2local(observer, c_set)
				moon_dist = int(np.max(sep_moon).value)
				priority = intbl['priority'][i]

				if header == False:
					f.write('#\tObservatory\t: {}\n'.format(observer.name))
					f.write('#\t{} UTC & Day Time Saving {}\n'.format('/'.join(time.value.split(' ')[0].split('-')), dts))
					f.write('#\t-18 deg sunset\t: {}\n'.format(ut2local(observer, sunset_tonight)))
					f.write('#\t-18 deg sunrise\t: {}\n'.format(ut2local(observer, sunrise_tonight)))
					# f.write('Moon ra, dec\t: \n')
					f.write('#\tMoon phase\t: {}\n'.format(moon_illumi))
					# f.write('#\tMoon seperation\t: {}\n'.format(np.max(sep_moon)))
					f.write('#\tMoon seperation limit\t>= {}\n'.format(constraint_moon))
					f.write('#\tAltitude limit\t>= {}\n'.format(constraint_alt))
					f.write('#{}\n'.format('-'*60))
					f.write('name ra dec rise(LT) transit(LT) set(LT) moon_dist(deg) priority\n')
					header = True
				line = '%-16s%-13s%-15s%-8s%-8s%-8s%-4s%-4s\n' % (name, ra, dec, rise_local, transit_local, set_local, moon_dist, priority)
				# print(line)
				f.write(line)
			else:
				# print('out of constraint', intbl[i]['obj'], 'moon.dist', np.max(sep_moon).value, 'alt', np.max(alt).value)
				pass
		except:
			# print('numpy.float64 item error', intbl[i]['obj'])
			pass
	f.close()
#============================================================
#	SETTING
#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_targetlist = '/home/sonic/Research/yourpy/gppy/table/imsng_alltarget.dat'
path_obs = '/home/sonic/Research/yourpy/gppy/table/observatory.txt'
path_out = '/home/sonic/Research/CEOU/IMSNG/targetlist/rts_vis_2021'
#------------------------------------------------------------
#	TABLE
obstbl = ascii.read(path_obs)
intbl = ascii.read(path_targetlist)
intbl = intbl[intbl['priority'] <= 2.0]
#------------------------------------------------------------
#	CONSTRAINT
constraint_alt = 30*u.deg		# [deg]
constraint_moon = 40*u.deg		# [deg]
#------------------------------------------------------------
#	OBSERVATORY INFO
obs = 'SAO'
# obs = 'LOAO'
# obs = 'McD'
observer = callobserver(obs, obstbl)
#------------------------------------------------------------
racol = intbl['ra']
decol = intbl['dec']
c = SkyCoord(racol, decol, unit=(u.hourangle, u.deg))

# obslist = ['SAO']
# obslist = ['LOAO']
# obslist = ['McD']

# y, m, d = 2020, 1, 1
y, m, d = 2021, 1, 1
# for obs in obslist:

for m_step in np.arange(0, 12, 1):
	m_input = m + m_step
	for d_step in np.arange(0, 32, 1):
		d_input = d + d_step
		try:
			param_rts = dict(	c = c,
								observer=observer,
								y=y, m=m_input, d=d_input,
								path_out=path_out)
			print('{}\t:{}/{}/{}'.format(obs, y, m_input, d_input))
			rts_vis_maker(**param_rts)
		except:
			print('{}/{}\tDONE'.format(y, m_input))


# datelist = [(2020, 1, 31), (2020, 2, 29), (2020, 3, 31), (2020, 4, 30), (2020, 5, 31), (2020, 6, 30), (2020, 7, 31), (2020, 8, 31), (2020, 9, 30), (2020, 10, 31), (2020, 11, 30), (2020, 12, 31), ]
datelist = [(2021, 1, 31), (2021, 2, 28), (2021, 3, 31), (2021, 4, 30), (2021, 5, 31), (2021, 6, 30), (2021, 7, 31), (2021, 8, 31), (2021, 9, 30), (2021, 10, 31), (2021, 11, 30), (2021, 12, 31), ]
for date in datelist:
	y, m, d = date[0], date[1], date[2]

	param_rts = dict(	c = c,
						observer=observer,
						y=y, m=m, d=d,
						path_out=path_out)
	print('{}\t:{}/{}/{}'.format(obs, y, m, d))
	rts_vis_maker(**param_rts)
