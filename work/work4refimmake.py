bandlist = []
for inim in photbl['obs']:
	bandlist.append(inim.split('-')[5])
photbl['band'] = np.array(bandlist)

photbl[ (photbl['band']=='R') & (photbl['seeing']<3.0) & (photbl['mag']>20.0)]
photbl[ (photbl['band']=='V') & (photbl['seeing']<4.0) & (photbl['mag']>19.0) ]
photbl[ (photbl['band']=='B') & (photbl['seeing']<2.5) & (photbl['mag']>20.6) ]

cp Calibrated-LOAO-NGC6814-20170702-080848-R-60.fits ./sel/
cp Calibrated-LOAO-NGC6814-20170702-081003-R-60.fits ./sel/
cp Calibrated-LOAO-NGC6814-20170702-081118-R-60.fits ./sel/
cp Calibrated-LOAO-NGC6814-20170702-091330-R-60.fits ./sel/
cp Calibrated-LOAO-NGC6814-20170702-091445-R-60.fits ./sel/
cp Calibrated-LOAO-NGC6814-20170702-091600-R-60.fits ./sel/



import os, glob

path_qso	= '/data3/IMSNG/LOAO/gal'
path_gundam	= '/mnt/window/Users/User/Downloads/data/loao/ref'
downcom		= 'sshpass -prjseka23! scp -ro StrictHostKeyChecking=no paek@qso.snu.ac.kr:'

#reflist		= ['NGC0488', 'NGC1309', 'UGC02855', 'NGC2207', 'NGC2993', 'IC2537', 'NGC3169', 'NGC3294', 'NGC3344', 'NGC3629', 'NGC3646', 'NGC3938', 'NGC4030', 'NGC4108']
reflist			= ['NGC3169', 'NGC3294', 'NGC3344', 'NGC3629', 'NGC3646', 'NGC3938']

for obj in reflist:
	os.system('mkdir '+path_gundam+'/'+obj+'/')
	com		= downcom+path_qso+'/'+obj+'/C*.fits '+path_gundam+'/'+obj+'/'
	print(com)
	os.system(com)