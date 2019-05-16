#============================================================
#	MAKE REGION FILE FOR ds9
#	20??.??.??	MADE BY		Changsu Choi
#	2019.05.16	UPDATED BY	Gregory S.H. Paek
#============================================================
import os,sys,glob
import string
from astropy.io import ascii 
import numpy as np
import math
#------------------------------------------------------------
def ds9reg(starname, ra, dec, outname='ds9', size=5, color='yellow'):
	radius	= """ {0}" """.format(size)
	color	= "yellow"
	os.system('rm '+outname+'.reg')
	f		= open(outname+'.reg','w')
	head1	= "# Region file format: DS9 version 4.1\n"; f.write(head1)
	head2	= """global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n"""; f.write(head2)
	head3	= "fk5\n"; f.write(head3)
	for n in range(len(name)):
		body="circle("+str(ra[n])+","+str(dec[n])+","+radius+") # color="+color+" text={"+str(starname[n])+"}\n"	
		f.write(body)
	f.close()
#------------------------------------------------------------
os.system('ls *.fits *.fits')
imlist	= glob.glob(input('IMAGE TO PROCESS\t: '))
