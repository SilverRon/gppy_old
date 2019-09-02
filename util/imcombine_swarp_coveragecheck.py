#	SIMPLE SWARP USAGE FOR COVERAGE CHECK
#	2019.08.30	CREATED BY Gregory S.H. Paek
#============================================================
#	MODULE
#------------------------------------------------------------
import os, sys, glob
#============================================================
obs = 'SAO'
inim_style = 'Calib-*-com.fits'
path_image = '/data1/S190425z/{}'.format(obs)
path_save = '/data1/S190425z/allcoverage'
imlist = glob.glob(path_image+'/'+inim_style)
imlistname = path_save+'/'+obs+'_swarplist.txt'
#------------------------------------------------------------
#	MAKE IMAGE LIST TO COMBINE
f = open(imlistname, 'a')
for inim in imlist:
	f.write(inim+'\n')
f.close()
#------------------------------------------------------------
os.sytem('swarp -dd > {}/default.swarp'.format(path_save))
os.system('swarp -c default.swarp @{}'.format(imlistname))
#------------------------------------------------------------
'''
import os, sys, glob
obs = 'KMTNet'
inim_style = 'a*.fits'
path_image = '/data1/S190425z/{}/raw'.format(obs)
path_save = '/data1/S190425z/allcoverage'
imlist = glob.glob(path_image+'/'+inim_style)
imlistname = path_save+'/'+obs+'_swarplist.txt'
f = open(imlistname, 'a')
for inim in imlist:
	f.write(inim+'\n')
f.close()
os.system('swarp -c default.swarp @{}'.format(imlistname))
'''


