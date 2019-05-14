#!/home/paek/anaconda2/bin/python2.7

#   1. FILE CHECK, whether new folder was updated.
#   2. pre-processing + log file
#   3. moving
#   4. combine and subtraction (not yet)
#   18.12.04 add something by G.Paek

import os, sys, glob
import numpy as np
import astropy.io.ascii as ascii

obs         = 'LOAO'

path_calib  = '/data3/paek/factory/calib/'
path_factory= '/data3/paek/factory/'
path_log    = '/home/paek/log/'
path_raw    = '/data3/grb/'
path_gal    = '/data3/IMSNG/IMSNGgalaxies/'
path_mframe = '/data3/paek/factory/master_frames/'

path_c_LOAO = '/data3/IMSNG/'+obs+'/2018/'
#------------------------------------------------------------------------------#
### 1. FOLDER CHECK
def folder_check(obs):                                # obs = 'loao'    
    '''
    CHECK OBSERVING DATE
    '''
    log     = obs.lower()+'.log'                    # loao.log
    log_read= ascii.read(path_log+log)              # /home/paek/log/loao.log
    date_lst= log_read['date']
    dir_lst = glob.glob(path_raw+obs.upper()+'/*/') # /data3/grb/LOAO/
    new     = []
    # IS THERE NEW FILE?
    for direc in dir_lst:
        if direc in date_lst: pass
        else:
            comment = 'New data\t: '+direc
            print(comment)
            new.append(direc)
    # MOVE NEW FILE TO FACTORY FOLDER
    for direc in new:
        directory   = direc.split('/')[-2]
        rmcom       = 'rm -rf '+path_factory+directory
        cpcom       = 'cp -r '+direc+' '+path_factory
        print(rmcom)    ; os.system(rmcom)
        print(cpcom)    ; os.system(cpcom)
    return new

#------------------------------------------------------------------------------#
import matplotlib.pyplot as plt
from astropy.modeling import models
from astropy import units as u
from astropy import nddata
from astropy.io import fits, ascii
from astropy.time import Time
from astropy.nddata import CCDData
import time
import astropy.coordinates as coord
from msumastro import ImageFileCollection, TableTree
import ccdproc

imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())

'''
CASE
C1		JD							SOAO, DOAO, BOAO, SAO1m, CCA250, T52(SNUCAMII)
C2		J_DATE						LOAO
C3		DATE-OBS					MAO(FLI), iTel(except T52)
C4		'DATE-OBS'+'T'+'TIME-OBS'	MAO(SNUCAM)			
C5		'DATE-OBS'+'T'+'UT'			30inch		
'''

###     DIRECTORY TO PROCESS
def call_images(path_folder):
    images          = ImageFileCollection(path_folder, keywords='*')
    return images

###     OBSERVATORY
obs         = 'LOAO'
path_info   = '/home/paek/qsopy/tables/obs.txt'
def obsinfo(obs, path_info):
    obs_table   = ascii.read(path_info)
    indx_obs    = np.where(obs == obs_table['obs'])
    gain        = float(obs_table['gain'][indx_obs][0]) * u.electron / u.adu
    rdnoise     = float(obs_table['RDnoise'][indx_obs][0]) * u.electron
    pixscale    = float(obs_table['pixelscale'][indx_obs][0])
    # dark        = float(obs_table['dark'][indx_obs][0])
    return gain, rdnoise, pixscale

#------------------------------------------------------------------------------#
### MASTER ZERO
def master_zero(images):
    start_time  = time.time()
    comment     = '-'*80+'\n' \
                + 'MAKING MASTER ZERO\n'
    print(comment)
    zero_list   = []
    for hdu, fname in images.hdus(imagetyp='Bias', return_fname=True):
        meta            = hdu.header
        meta['FILENAME']= fname
        zero_list.append(ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu"))
    comment     = str(len(zero_list))+' zero images will be combined.'  ; print(comment)
    zeros       = ccdproc.Combiner(zero_list)
    mzero       = zeros.median_combine()
    # header
    newmeta     = meta
    for n in range(0, len(zero_list)):
        newkey  = 'COMB'+str(n+1)
        newvalue= zero_list[n].header['FILENAME']
        newmeta[newkey] = newvalue
    newmeta['FILENAME'] = 'zero.fits'
    mzero.header  = newmeta
    os.system('rm '+new+'/zero.fits')
    mzero.write(new+'/zero.fits')
    zero_min, zero_max, zero_mean, zero_std = imstats(np.asarray(mzero))
    plt.figure(figsize=(15, 15))
    plt.imshow(mzero, vmax=zero_mean + 4*zero_std, vmin=zero_mean - 4*zero_std)
    plt.savefig(new+'/zero.png')
    times   = round(time.time() - start_time)
    print(str(times)+' seconds')
    return  mzero
#------------------------------------------------------------------------------#
### MASTER FLAT
def master_flat(images, mzero, filte):
    start_time  = time.time()
    comment     = '-'*80+'\n' \
                + 'MAKING MASTER FLAT\n'
    print(comment)
    comment = 'Subtract zero from '+filte+' flat images.'    ; print(comment)
    outname     = 'n'+filte+'.fits'
    zflat_list  = []
    for hdu, fname in images.hdus(imagetyp='flat', filter=filte, return_fname=True):
        meta    = hdu.header
        meta['filename'] = fname
        flat    = ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu")
        zflat   = ccdproc.subtract_bias(flat, mzero)
        subhead = 'ccd='+fname+', master='+mzero.meta['FILENAME']
        zflat.meta['SUBBIAS']   =subhead
        zflat_list.append(ccdproc.CCDData(data=zflat.data, meta=meta, unit="adu"))
    comment     = str(len(zflat_list))+' '+filte+' flat images will be combined.'  ; print(comment)
    flat_combiner   = ccdproc.Combiner(zflat_list)
    flat_combiner.minmax_clipping()
    scaling_func    = lambda arr: 1/np.ma.median(arr)
    flat_combiner.scaling   = scaling_func

    mflat     = flat_combiner.median_combine()
    newmeta         = meta
    for n in range(0, len(zflat_list)):
        newkey  = 'COMB'+str(n+1)
        newvalue= zflat_list[n].header['FILENAME']
        newmeta[newkey] = newvalue
    newmeta['FILENAME'] = outname
    mflat.header  = newmeta
    os.system('rm '+new+'/'+outname)
    mflat.write(new+'/'+outname)
    f_min, f_max, f_mean, f_std = imstats(np.asarray(mflat))
    plt.figure(figsize=(15, 15))
    plt.imshow(mflat, vmin=f_mean-5*f_std, vmax=f_mean+5*f_std)
    plt.savefig(new+'/'+outname[:-5]+'.png')
    mflat_electron = ccdproc.gain_correct(mflat, gain=gain)
    times   = round(time.time() - start_time)
    print(str(times)+' seconds')
    return mflat
#------------------------------------------------------------------------------#
### OBJECT - MASTER BIAS
def calibration(images, mzero, mflat, filte):
    start_time  = time.time()
    comment     = '-'*80+'\n' \
                + 'OBJECT CORRECTION\n'
    print(comment)
    for hdu, fname in images.hdus(imagetyp='Light', filter=filte, return_fname=True):
        meta            = hdu.header
        meta['filename']= fname
        obj_list.append(meta['object'])
        obj     = ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu")
        # ZERO CORRECTION
        zobj    = ccdproc.subtract_bias(obj, mzero)
        subhead = 'ccd='+fname+', master='+mzero.meta['FILENAME']
        meta['SUBBIAS'] = subhead
        # FLAT CORRECTION
        fzobj   = ccdproc.flat_correct(zobj, mflat)
        divhead = 'ccd='+fname+', master='+mflat.meta['FILENAME']
        meta['DIVFLAT'] = divhead
        fzobj.header    = meta
        fzobj.write(new+'/'+'fz'+fname)
    times   = round(time.time() - start_time)
    print(str(times)+' seconds')

    
#------------------------------------------------------------------------------#
### ASTROMETRY
def astrometry(direc, pixscale):
    start_time  = time.time()
    order   = 0
    imlist  = glob.glob(direc+'/fz*.fits')
    upscl   = str(pixscale + pixscale*0.05)
    loscl   = str(pixscale - pixscale*0.05)
    comment = '='*80+'\n' \
            + 'ASTROMETRY PROCESS START\t'+str(len(imlist))+' IMAGES'
    print(comment)
    for im  in imlist:
        order   += 1
        img     = im.split('/')[-1]
        outname = new+'/a'+img
        com     ='solve-field '+im \
                +' --scale-unit arcsecperpix --scale-low '+loscl+' --scale-high '+upscl \
                +' --no-plots --new-fits '+outname+' --overwrite --use-sextractor\n'
        os.system(com)
        comment = '-'*80+'\n' \
                + 'ASTROMETRY PROCESS\t'+str(round(order*100/len(imlist)))+'%\t'\
                + '['+str(order)+'/'+str(len(imlist))+']'+'\n'
        print(comment)
    os.system('rm '+direc+'/*.axy '+direc+'/*.corr '+direc+'/*.xyls '+direc+'/*.match '+direc+'/*.rdls '+direc+'/*.solved '+direc+'/*.wcs ')
    comment     = 'COMPLETE\n'+'='*80
    print(comment)
    times   = round(time.time() - start_time)
    print(str(times)+' seconds')
#------------------------------------------------------------------------------#
### PUT MJD IN HEADER
def putMJD(inim, case):
	hdr=fits.getheader(inim)
	if   case == 'C1' : 
		mjdval=Time(hdr['JD'],format='jd',scale='utc').mjd
		fits.setval(inim,'MJD',value=mjdval,comment='MJD appended')	
	elif case == 'C2' : 
		mjdval=Time(float(hdr['J_DATE'][1:-24]),format='jd',scale='utc').mjd
		fits.setval(inim,'MJD',value=mjdval,comment='MJD appended')	
	elif case == 'C3' :
		mjdval=Time(hdr['DATE-OBS'],format='isot',scale='utc').mjd
		fits.setval(inim,'MJD',value=mjdval,comment='MJD appended')	
	elif case == 'C4' : 
		datetimestr=hdr['DATE-OBS']+'T'+hdr['TIME-OBS']
		mjdval = Time(datetimestr,format='isot',scale='utc').mjd
		fits.setval(inim,'MJD',value=mjdval,comment='MJD appended')	
	elif case =='C5' : 
		datetimestr=hdr['DATE-OBS']+'T'+hdr['UT']
		mjdval = Time(datetimestr,format='isot',scale='utc').mjd
		fits.setval(inim,'MJD',value=mjdval,comment='MJD appended')	
	else : print 'Not in known cases.'
	comment     = inim+'\t'+': '+'Header Updated to '+str(fits.getheader(inim)['MJD'])+'\t'+'[MJD]'
	#print(comment)
#------------------------------------------------------------------------------#
def isot_to_mjd(time):      # 20181026 to 2018-10-26T00:00:00:000 to MJD form
    yr  = time[0:4]         # year
    mo  = time[4:6]         # month
    da  = time[6:8]         # day
    isot = yr+'-'+mo+'-'+da+'T00:00:00.000'         # ignore hour:min:sec
    t = Time(isot, format='isot', scale='utc').mjd  # transform to MJD
    return t
#------------------------------------------------------------------------------#
def fnamechange(inim, obs) :
    hdr     = fits.getheader(inim)
    dateobs = hdr['DATE-OBS']
    datestr = dateobs[0:4]+dateobs[5:7]+dateobs[8:10]
    timestr = dateobs[11:13]+dateobs[14:16]+dateobs[17:19]
    objname = hdr['OBJECT']
    objname = objname.upper()	
    fil     = hdr['FILTER']  # str(hdr['FILTER'])
    if  fil     == 'B102'   : filname = 'B'
    elif fil    == 'V103'   : filname = 'V'
    elif fil    == 'R104'   : filname = 'R'
    elif fil    == 'I105'   : filname = 'I'
    elif fil    == 'Z'      : filname = 'Z'
    elif fil    == 'Y'      : filname = 'Y'	
    else                    : filname = 'NA'
    exptime = str(int(hdr['EXPTIME']))
    newname = 'Calib-'+obs+'-'+objname+'-'+datestr+'-'+timestr+'-'+filname+'-'+exptime+'.fits'
    cpcom   = 'cp '+inim+' '+new+'/'+newname
    os.system(cpcom)
    comment = newname+' is generated.'
    print(comment)
#------------------------------------------------------------------------------#
def move_calib(obs, path_calib, path_gal):
    list_gal    = []
    for img in glob.glob(path_calib+'Calib-'+obs+'-*.fits'):
        tmp     = img.split('-')[2]
        list_gal.append(tmp)
    list_gal    = list(set(list_gal))
    ##  LIST FOR OBJECT LIST
    list_dir    = []
    for obj in glob.glob(path_gal+'*/'):
        tmp     = obj.split('/')[-2]
        list_dir.append(tmp)
    list_dir    = list(set(list_dir))
    ##  MOVE FILES
    list_moved  = []
    list_created= []
    for gal in list_gal:
        if gal in list_dir:
            list_moved.append(gal)
            com = 'mkdir '+path_gal+gal+'/'+obs+'/'
            os.system(com)
            print(com)
            com = 'cp '+path_calib+'*'+gal+'*.fits '+path_gal+gal+'/'+obs+'/'
            os.system(com)
            print(com)
        else:
            list_created.append(gal)
            com = 'mkdir '+path_gal+gal+'/'
            os.system(com)
            print(com)
            com = 'mkdir '+path_gal+gal+'/'+obs+'/'
            os.system(com)  ; print(com)
            com = 'cp '+path_calib+'*'+gal+'*.fits '+path_gal+gal+'/'+obs+'/'
            os.system(com)  ; print(com)
        rmcom   = 'rm '+path_calib+'*'+gal+'*.fits'
        print(rmcom)  ; os.system(rmcom)
    comment     = '-'*80
    print(comment)
    comment     = '\nBELOW OBJECTS IMAGES WERE MOVED.'  ; print(comment)
    for obj in list_moved: print(obj)
    comment     = '\nBELOW OBJECTS FOLDER WERE GENERATED AND MOVED.'    ; print(comment)
    for obj in list_created: print(obj)
#==============================================================================#
### CALIBRATION
obs         = 'LOAO'
case        = 'C2'
path_info   = '/home/paek/qsopy/tables/obs.txt'
gain, rdnoise, pixscale = obsinfo(obs, path_info)
new_LOAO    = folder_check(obs)
### CALIBRATION
order       = 0
for path in new_LOAO:
    order   += 1
    new     = path_factory+path.split('/')[-2]
    date    = new.split('/')[-1]    ; date  = ('').join(date.split('_'))
    comment = 'REDUCTION START\t['+str(order)+'/'+str(len(new_LOAO))+'] : '+path; print(comment)
    images  = call_images(new)
    # MASTER ZERO
    mzero   = master_zero(images)
    cpcom   = 'cp '+new+'/zero.fits '+path_mframe+obs+'/zero/'+date+'-zero.fits'
    print(cpcom)    ; os.system(cpcom)
    cpcom   = 'cp '+new+'/zero.png '+path_mframe+obs+'/zero/'+date+'-zero.png'
    print(cpcom)    ; os.system(cpcom)

    #------------------------------------------------------------------------------#
    # OBJECT & FILTER LIST OF SCI IMAGES
    obj_list= []
    fil_list= []
    for hdu, fname in images.hdus(imagetyp='Light', return_fname=True):
        meta= hdu.header
        fil_list.append(meta['FILTER'])
        obj_list.append(meta['OBJECT'])
    obj_list= list(set(obj_list))
    fil_list= list(set(fil_list))
    #------------------------------------------------------------------------------#
    # FILTER LIST OF FLAT IMAGES
    band_list= []
    for hdu, fname in images.hdus(imagetyp='Flat', return_fname=True):
        meta= hdu.header
        band_list.append(meta['FILTER'])
    band_list= list(set(band_list))
    
    order_fil   = 0
    for filte in fil_list:
        # WHETHER FLAT IS ENOUGHT
        if filte not in band_list:
            # IF THERE IS NO FLAT FRAMES, BORROW FROM CLOSEST OTHER DATE
            order_fil   += 1
            comment = '\n' + '-'*80 + '\n' \
                    + 'No flat frames for '+filte + '\n' \
                    + '-'*80 + '\n'
            print(comment)
            mframe_flat = np.array( glob.glob(path_mframe+obs+'/flat/*n'+filte+'.fits') )
            mframe_zero = np.array( glob.glob(path_mframe+obs+'/zero/*zero.fits') )
            diff        = []
            for flat in mframe_flat:
                onlyfile    = flat.split('/')[-1]
                dateinfo    = onlyfile.split('-')[0]
                diff.append(isot_to_mjd(date)-isot_to_mjd(dateinfo))
            indx_closet = np.where(diff == np.min(diff))
            cpcom   = 'cp '+mframe_flat[indx_closet][0]+' '+new+'/'
            print(cpcom)    ; os.system(cpcom)
            mflat   = CCDData.read(mframe_flat[indx_closet][0], hdu=0, unit='adu')
        #------------------------------------------------------------------------------#
        else:
            # IF THERE IS A FLAT FRAMES, JUST DO IT
            # MASTER FLAT
            order_fil   += 1
            comment = 'Processing\t'+filte+'\t['+str(order)+'/'+str(len(fil_list))+']'  ; print(comment)
            mflat   = master_flat(images, mzero, filte)
            cpcom   = 'cp '+new+'/n'+filte+'.fits '+path_mframe+obs+'/flat/'+date+'-n'+filte+'.fits'
            print(cpcom)    ; os.system(cpcom)
            cpcom   = 'cp '+new+'/n'+filte+'.png '+path_mframe+obs+'/flat/'+date+'-n'+filte+'.png'
            print(cpcom)    ; os.system(cpcom)
            #------------------------------------------------------------------------------#
        calibration(images, mzero, mflat, filte)
    #------------------------------------------------------------------------------#
    astrometry(new, pixscale)

    imlist  = glob.glob(new+'/afz*.fits')
    for img  in imlist:
#        img     = im.split('/')[-1]
        putMJD(img, case)
        fnamechange(img, obs)
    os.system('rm '+new+'/*fz*.fits')
    os.system('chmod 777 '+new+'/Calib*.fits')
    os.system('cp '+new+'/Calib*.fits '+path_calib)
    os.system('rm '+new+'/Calib*.fits')
    os.system('cp -rf '+new+'/ '+path_c_LOAO)
    os.system('rm -rf '+new+'/')
    f       = open(path_log+'loao.log', 'a')
    f.write(path+'\n')
    f.close()
    comment = 'REDUCTION DONE\t['+str(order)+'/'+str(len(new_LOAO))+'] : '+path+'\n' \
            + '='*80+'\n'
    print(comment)
#------------------------------------------------------------------------------#
### REMOVE UNNESSARY FILES AND MOVE CALIBRATED FILES
##  LIST FOR CALIB IMAGES
if len(new_LOAO) != 0: move_calib(obs, path_calib, path_gal)

