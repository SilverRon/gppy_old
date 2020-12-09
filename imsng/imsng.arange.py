import os
import glob

#   Arangement images (format : Calib-[obs]-[object]-[20XXXXXX]-[XXXXXX]-[filter]-[exptime].fits)

imlist = sorted(glob.glob('Calib*.fits'))
filterlist = []

for inim in imlist:
    part = inim.split('-')
    filte = part[-2]
    filterlist.append(filte)

for filte in set(filterlist):
    com0 = 'mkdir {}'.format(filte)
    com1 = 'mv ./Calib*-{}-*.fits ./{}/'.format(filte, filte)
    os.system(com0)
    print(com0)
    os.system(com1)
    print(com1)

os.system('ls ./*/ -d')
