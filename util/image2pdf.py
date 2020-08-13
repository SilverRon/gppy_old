#	Image file to pdf file FOR PYTHON 3
#	CREATED	2020.08.07	Gregory S.H. Paek
#	UPDATE 
#============================================================
def img2pdf(inim, outim):
	'''
	inim = './test.png'
	outim = './test.pdf'
	'''
	from PIL import Image
	image1 = Image.open(inim)
	im1 = image1.convert('RGB')
	im1.save(outim)