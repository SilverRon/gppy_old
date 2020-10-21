#%%
#	MULTI PLOT SAMPLE
#============================================================
import glob
import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.time import Time
import matplotlib.pyplot as plt
from ligo.skymap.plot.marker import reticle
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
from astropy.visualization import ZScaleInterval, LinearStretch
from matplotlib.patches import Circle
from imsng import tool

path_save = '.'

plt.rc('font', family='serif')
plt.close('all')
fig = plt.figure()

ax0 = fig.add_subplot(111)
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
# ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)
# ax6 = fig.add_subplot(236)

x = 1920 / fig.dpi
y = 1080 / fig.dpi
fig.set_figwidth(x)
fig.set_figheight(y)


ax0.spines['top'].set_color('none')
ax0.spines['bottom'].set_color('none')
ax0.spines['left'].set_color('none')
ax0.spines['right'].set_color('none')
ax0.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

ax0.plot(1, 1, color='tomato', label='Test1')
ax0.plot(1, 1, color='dodgerblue', label='Test2')
ax0.plot(1, 1, color='gold', label='Test3')

ax0.legend(fontsize=14, ncol=3)

fig.savefig("{}/multi-plot.pdf".format(path_save), bbox_inches='tight', overwrite=True)
fig.savefig("{}/multi-plot.png".format(path_save), bbox_inches='tight', dpi=500, overwrite=True)