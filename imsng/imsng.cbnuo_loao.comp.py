#	This code is for drawing figure in pdf extinction for Paper
#============================================================
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.time import Time
#============================================================
#	Function
#------------------------------------------------------------
def imsng_dist(m):
	d = 50.0 * (10.**((m-19.5)/5.))
	return d
#============================================================
#	User setting
#------------------------------------------------------------
#	PATH
#------------------------------------------------------------
path_loao = '/data1/IMSNG/loao/2020_1017/phot.dat'
path_cbnuo = '/data1/IMSNG/cbnuo/20201018-imsng/phot.dat'

path_loao_com_R = '/data1/IMSNG/loao/2020_1017/com/R/phot.dat'
path_loao_com_B = '/data1/IMSNG/loao/2020_1017/com/B/phot.dat'
path_cbnuo_com = '/data1/IMSNG/cbnuo/20201018-imsng/com/R/phot.dat'
path_doao_com = '/data1/IMSNG/doao/20201017-1m-IMSNG/com/phot.dat'

path_save = '/home/sonic/Pictures'
#------------------------------------------------------------
ltbl = ascii.read(path_loao)
loaotbl = ltbl
loaotbl = loaotbl[np.argsort(loaotbl['jd'])]
ltbl = ltbl[np.argsort(ltbl['jd'])]
ltbl = ltbl[ltbl['filter']=='R']


ctbl = ascii.read(path_cbnuo)
ctbl = ctbl[np.argsort(ctbl['jd'])]
lcmtbl = ascii.read(path_loao_com_R)
lcmtbl = lcmtbl[np.argsort(lcmtbl['jd'])]

lcmtbl_ = ascii.read(path_loao_com_B)
lcmtbl_ = lcmtbl_[np.argsort(lcmtbl_['jd'])]

lcmtbl_all = vstack([lcmtbl, lcmtbl_])

ccmtbl = ascii.read(path_cbnuo_com)
ccmtbl = ccmtbl[np.argsort(ccmtbl['jd'])]

dcmtbl = ascii.read(path_doao_com)
dcmtbl = dcmtbl[np.argsort(dcmtbl['jd'])]

dcmtbl_ = dcmtbl[dcmtbl['filter']=='R']
#============================================================
#	PLOT
#------------------------------------------------------------
plt.rc('font', family='serif')
plt.close('all')
fig = plt.figure()
ax0 = fig.add_subplot(111)
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(212)
# ax4 = fig.add_subplot(224)
#------------------------------------------------------------
x = 1920 / fig.dpi
y = 1080 / fig.dpi
fig.set_figwidth(x)
fig.set_figheight(y)
#------------------------------------------------------------
ax0.spines['top'].set_color('none')
ax0.spines['bottom'].set_color('none')
ax0.spines['left'].set_color('none')
ax0.spines['right'].set_color('none')
ax0.tick_params(labelcolor='w', top=False,
                bottom=False, left=False, right=False)
ax0.set_xlabel('JD', fontsize=20)
# ax0.set_ylabel('Seeing [arcsec]', fontsize=20)
#============================================================
#	SEEING
#------------------------------------------------------------
#	LOAO
#------------------------------------------------------------
xl = np.arange(np.min(ltbl['jd']), np.max(ltbl['jd']), 0.01)
yl_seeing = np.median(ltbl['seeing'])
yl_seeing_up = yl_seeing + np.std(ltbl['seeing'])
yl_seeing_lo = yl_seeing - np.std(ltbl['seeing'])
#------------------------------------------------------------
ax1.plot(ltbl['jd'], ltbl['seeing'], marker='o', linestyle='-', color='gold', alpha=0.5)
# ax1.fill_between(xl, yl_seeing_lo, yl_seeing_up, color='gold', alpha=0.5)
ax1.axhline(y=yl_seeing, color='k', linestyle='--', label='Median : {}'.format(round(yl_seeing, 1)))
ax1.set_ylabel('Seeing [arcsec]', fontsize=20)
ax1.legend(fontsize=20)
#------------------------------------------------------------
#	CBNUO
#------------------------------------------------------------
xc = np.arange(np.min(ctbl['jd']), np.max(ctbl['jd']), 0.01)
yc_seeing = np.median(ctbl['seeing'])
yc_seeing_up = yc_seeing + np.std(ctbl['seeing'])
yc_seeing_lo = yc_seeing - np.std(ctbl['seeing'])
#------------------------------------------------------------
ax2.plot(ctbl['jd'], ctbl['seeing'], marker='o', linestyle='-', color='green', alpha=0.5)
# ax2.fill_between(xc, yl_seeing_lo, yl_seeing_up, color='gold', alpha=0.5)
ax2.axhline(y=yc_seeing, color='k', linestyle='-.', label='Median : {}'.format(round(yc_seeing, 1)))
ax2.legend(fontsize=20)

ax1.tick_params(axis="x", labelsize=14)
ax1.tick_params(axis="y", labelsize=14)
ax2.tick_params(axis="x", labelsize=14)
ax2.tick_params(axis="y", labelsize=14)

#============================================================
#	DEPTH
#------------------------------------------------------------
#	LOAO
#------------------------------------------------------------
xl = np.arange(np.min(ltbl['jd']), np.max(ltbl['jd']), 0.01)
yl_depth = np.median(ltbl['ul_5sig'])
yl_depth_up = yl_depth + np.std(ltbl['ul_5sig'])
yl_depth_lo = yl_depth - np.std(ltbl['ul_5sig'])
#------------------------------------------------------------
ax3.plot(ltbl['jd'], ltbl['ul_5sig'], marker='o', linestyle='-', color='gold', alpha=0.5)
ax3.set_ylabel(r'5$\sigma$ depth [AB mag]', fontsize=20)
#------------------------------------------------------------
#	COM
#------------------------------------------------------------
ax3.axhline(y=np.median(lcmtbl['ul_5sig']), color='darkorange', linestyle='--', linewidth=2.5, label='LOAO : {}'.format(round(np.median(lcmtbl['ul_5sig']), 1)))
ax3.plot(lcmtbl['jd'], lcmtbl['ul_5sig'], marker='v', ms=10, mec='k', linestyle='', color='gold', alpha=1.0, label='60s*3')
#------------------------------------------------------------
#	CBNUO
#------------------------------------------------------------
xc = np.arange(np.min(ctbl['jd']), np.max(ctbl['jd']), 0.01)
yc_depth = np.median(ctbl['ul_5sig'])
yc_depth_up = yc_depth + np.std(ctbl['ul_5sig'])
yc_depth_lo = yc_depth - np.std(ctbl['ul_5sig'])
#------------------------------------------------------------
ax3.plot(ctbl['jd'], ctbl['ul_5sig'], marker='o', linestyle='-', color='green', alpha=0.5)
#------------------------------------------------------------
#	COM
#------------------------------------------------------------
ax3.axhline(y=np.median(ccmtbl['ul_5sig']), color='green', linestyle='-.', linewidth=2.5, label='CBNUO : {}'.format(round(np.median(ccmtbl['ul_5sig']), 1)))
ax3.plot(ccmtbl['jd'], ccmtbl['ul_5sig'], marker='v', ms=10, mec='k', linestyle='', color='green', alpha=1.0, label='180s*5')
#------------------------------------------------------------
#	DOAO
#------------------------------------------------------------
ax3.axhline(y=np.median(dcmtbl_['ul_5sig']), color='blue', linestyle=':', linewidth=2.5, label='DOAO : {}'.format(round(np.median(dcmtbl['ul_5sig']), 1)))
#------------------------------------------------------------
d_loao, u_loao = ax3.set_ylim()
d_cbnuo, u_cbnuo = ax3.set_ylim()
ax3.set_ylim([u_loao, d_cbnuo])
ax3.legend(fontsize=18)
# ax4.set_ylim([up, down])
# ax4.legend(fontsize=20)

ax3.tick_params(axis="x", labelsize=14)
ax3.tick_params(axis="y", labelsize=14)
#------------------------------------------------------------
ax33 = ax3.twinx()
ax33.set_ylabel('D [Mpc]', fontsize=20)
ax33.set_ylim([imsng_dist(u_loao), imsng_dist(d_cbnuo)])
ax33.set_yscale('log')
yticks = np.arange(10, 60, 10)
ax33.set_yticks(yticks)
ax33.set_yticklabels(yticks)
ax33.tick_params(axis="y", labelsize=14)
ax33.minorticks_on()
#============================================================
#	SAVE
#------------------------------------------------------------
fig.savefig("{}/imsng.cbnuo_loao.comp.pdf".format(path_save), bbox_inches='tight', overwrite=True)
fig.savefig("{}/imsng.cbnuo_loao.comp.png".format(path_save), bbox_inches='tight', overwrite=True)
#============================================================
#	
# lcmtbl_all = lcmtbl_all[np.argsort(lcmtbl_all['jd'])]
# lcmtbl_all['t-t0'] = lcmtbl_all['jd']-np.min(lcmtbl_all['jd'])
# delt_loao = lcmtbl_all['t-t0'][1:] - lcmtbl_all['t-t0'][:-1]
lcmtbl = lcmtbl[np.argsort(lcmtbl['jd'])]
lcmtbl['t-t0'] = lcmtbl['jd']-np.min(lcmtbl['jd'])
delt_loao = lcmtbl['t-t0'][1:] - lcmtbl['t-t0'][:-1] #-3/(24*60)
#	[min]
# print('LOAO')
# print(np.mean(delt_loao)*24*60, np.std(delt_loao)*24*60)
# print(np.min(delt_loao)*24*60, np.max(delt_loao)*24*60)
ccmtbl['t-t0'] = ccmtbl['jd']-np.min(ccmtbl['jd'])
delt_cbnuo = ccmtbl['t-t0'][1:] - ccmtbl['t-t0'][:-1]
#	[min]
# print('CBNUO')
# print(np.mean(delt_cbnuo)*24*60, np.std(delt_cbnuo)*24*60)
# print(np.min(delt_cbnuo)*24*60, np.max(delt_cbnuo)*24*60)
dcmtbl_['t-t0'] = dcmtbl_['jd']-np.min(dcmtbl_['jd'])
delt_doao = dcmtbl_['t-t0'][1:] - dcmtbl_['t-t0'][:-1]
#------------------------------------------------------------





bins = np.arange(0, 82.5, 2.5)

plt.close('all')
plt.hist(delt_loao*24*60, bins=bins, color='gold', alpha=0.5, label='LOAO BR ({}) {} min'.format(int(len(lcmtbl_all)/2), round(np.mean(delt_loao)*24*60, 1)))
plt.hist(delt_doao*24*60, bins=bins, color='dodgerblue', alpha=0.5, label='DOAO BVR ({}) {} min'.format(int(len(dcmtbl)), round(np.mean(delt_doao)*24*60, 1)))
plt.hist(delt_cbnuo*24*60, bins=bins, color='green', alpha=0.5, label='CBNUO R ({}) {} min'.format(int(len(ccmtbl)), round(np.mean(delt_cbnuo)*24*60, 1)))


params_plot = dict(	
					x=[np.mean(delt_cbnuo)*24*60],
					y=[15],
					xerr=[np.std(delt_cbnuo)*24*60],
					marker='s', ms=15, mew=1,
					mec='k', mfc='green', alpha=1.0,
					fmt='ko', 
					capsize=6.5, capthick=2,
				    # label='CBNUO {} min'.format(round(np.mean(delt_cbnuo)*24*60, 0))
					)
plt.errorbar(**params_plot)

params_plot = dict(
					x=[np.mean(delt_loao)*24*60],
					y=[15],
					xerr=[np.std(delt_loao)*24*60],
					marker='s', ms=15, mew=1,
					mec='k', mfc='gold', alpha=1.0,
								fmt='ko',
					capsize=6.5, capthick=2,
					# label='CBNUO {} min'.format(round(np.mean(delt_cbnuo)*24*60, 0))
				)
plt.errorbar(**params_plot)

params_plot = dict(
					x=[np.mean(delt_doao)*24*60],
					y=[10],
					xerr=[np.std(delt_doao)*24*60],
					marker='s', ms=15, mew=1,
					mec='k', mfc='dodgerblue', alpha=1.0,
								fmt='ko',
					capsize=6.5, capthick=2,
					# label='CBNUO {} min'.format(round(np.mean(delt_cbnuo)*24*60, 0))
				)
plt.errorbar(**params_plot)


plt.legend(fontsize=18, loc='upper right', framealpha=1.0)
plt.xlabel('observation time per target [min]', fontsize=20)
plt.ylabel('#', fontsize=20)
# plt.xscale('log')
left, right = plt.xlim()
down, up  = plt.ylim()
plt.xlim(0, right)
plt.ylim(0, up)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid('both', linestyle='--', color='grey', alpha=0.5)
plt.tight_layout()


plt.savefig("{}/imsng.overhead.pdf".format(path_save), bbox_inches='tight', overwrite=True)
plt.savefig("{}/imsng.overhead.png".format(path_save), dpi=500, bbox_inches='tight', overwrite=True)
