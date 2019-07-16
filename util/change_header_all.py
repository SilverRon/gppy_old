import numpy as np
from astropy.io import fits
import os,sys,glob

# file_name = raw_input('FILE NAME : ')
list_name = raw_input('IMAGE TO PROCESS\t: ')
where = raw_input('WHICH HEADER\t: ')
what = raw_input('TO WHAT\t: ')


# input file list
# imglist = np.genfromtxt(list_name, usecols=(0), dtype=str)
# imglist = list(imglist)
imglist = glob.glob(list_name)

#for multiple data
def chg_header(imglist, where, what) :
	for img in imglist :
		data, header = fits.getdata(img, header = True)
		fits.getheader(img, 0)
		header[where] = what
		fits.writeto(img, data, header, clobber=True)

chg_header(imglist, where, what)

for i in imglist : os.system('gethead '+i+' '+where)

'''
data, header = fits.getdata(file_name, header = True)
fits.getheader(file_name, 0)
header[where] = what
fits.writeto(file_name, data, header, clobber=True)
'''
print 'DONE'
