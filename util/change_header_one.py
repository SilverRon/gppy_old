file_name = raw_input('FILE NAME : ')
where = raw_input('WHICH HEADER : ')
what = raw_input('TO WHAT : ')

#for multiple data
def chg_header(imglist, where, what) :
	for img in imglist :
		data, header = fits.getdata(img, header = True)
		fits.getheader(img, 0)
		header[where] = what
		fits.writeto(img, data, header, clobber=True)

from astropy.io import fits
data, header = fits.getdata(file_name, header = True)

fits.getheader(file_name, 0)
header[where] = what
fits.writeto(file_name, data, header, clobber=True)