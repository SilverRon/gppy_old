import os, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import warnings
warnings.filterwarnings("ignore")
#============================================================
#	FUNCTION
#------------------------------------------------------------
def img2mp4(paths, pathOut, fps=10):
	import cv2
	frame_array = []
	for idx, path in enumerate(paths):
		img = cv2.imread(path)
		height, width, layers = img.shape
		size = (width, height)
		frame_array.append(img)
	out = cv2.VideoWriter(pathOut, cv2.VideoWriter_fourcc(*'DIVX'), fps, size)
	for i in range(len(frame_array)):
		# writing to a image array
		out.write(frame_array[i])
	out.release()


#============================================================
#	EXAMPLE
#------------------------------------------------------------
# paths = sorted(glob.glob('/data1/S190425z/200625.result/skymap/*.png'))
# img2mp4(paths, './test.mp4', fps=10)
