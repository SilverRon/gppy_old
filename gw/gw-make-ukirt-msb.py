from astropy.io import ascii
import numpy as np
'''
PKS0106         01 08 38.771 + 01 35 00.32 RJ
3C84            03 19 48.160 + 41 30 42.10 RJ
3C120           04 33 11.096 + 05 21 15.62 RJ
PKS0438         04 40 17.180 - 43 33 08.60 RJ
'''

intbl = ascii.read('/mnt/window/Users/User/Downloads/data/Project/gw/S190426c/S190426c_Update-20190426-rts_vis-UKIRT-forOT.txt')

uname, ura, udec = intbl['name'], intbl['ra'], intbl['dec'] 

for i in range(0, 50):
	rapart = ura[i].split(':')
	depart = udec[i].split(':')

	newra = ' '.join(rapart)
	
	if float(depart[0]) > 0:
		depart.insert(0, '+')
	elif float(depart[0]) < 0:
		depart[0] = str(abs(int(depart[0])))
		depart.insert(0, '-')

	newdec = ' '.join(depart)

	line = uname[i]+' '+newra+' '+newdec+' RJ'
	print(line)

