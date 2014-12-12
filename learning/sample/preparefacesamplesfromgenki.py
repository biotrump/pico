#python3
# 1. reading the image file list and covert the color RGB image file to a 8bit grey scale image.
# 2. cropping a new box around the face bounding boxes. The new cropping area is no bigger than 192*192.
# 3. save the new cropping box as a rid file.
# 4. rid file : 4 bytes width, 4 bytes height, row-wise binary pixel data
#
import random
import numpy
import matplotlib.pyplot
import matplotlib.image
import matplotlib.cm
#import Image #python2
#python3
from PIL import Image
from PIL import ImageOps
import struct
import argparse
import os

#parsing command line args to src folder and dst folder
parser = argparse.ArgumentParser()
parser.add_argument('srcfolder', help='input folder, source')
parser.add_argument('dstfolder', help='output folder, destination')
args = parser.parse_args()
print ('args:',args)

#get the src and target folder from command line args
srcfolder = args.srcfolder
dstfolder = args.dstfolder
print ('src folder:',srcfolder)
print ('dst folder:',dstfolder)

# create destination folder, if needed
if not os.path.exists(dstfolder):
   os.makedirs(dstfolder)

#debug to enable plot the picture with a face bounding box
plot = 0

if plot:
	fig = matplotlib.pyplot.figure()
	matplotlib.pyplot.show(block=False)

#exit()
#im is an array/list from a cropped image around the face bounding box:
#im dim : row * col
#rid file : 4 bytes width, 4 bytes height , image data (row-wise)
def saveasrid(im, path):
	#
	# raw intensity data
	#

	#im is a array/list type: row(height) * col(width)
	h = im.shape[0]
	w = im.shape[1]
	print('h=',h,' w=',w)
	#
	f = open(path, 'wb')

	#The first record of rid file is binary w,h
	#'ii' = 2 (i)nteger (w,h) to write, each interger is 4 bytes, LSB
	data = struct.pack('ii', w, h)
	f.write(data)
	#create an array/list w*h with all the value "None" which means "null", "empty", False, "0"
	tmp = [None]*w*h #[None,None,None.......]
	#converting 2D of im to 1D tmp
	for y in range(0, h):
		for x in range(0, w):
			tmp[y*w + x] = im[y, x]

	#write binary data of im
	#http://www.tutorialspoint.com/python/python_strings.htm
	#'%sB' % w*h : refer to format string, %s is used for w*h to convert to string type before B.
	#for example: if w*h= 5*6=30, '%sB' % w*h ='%sB' % 5*6 = '%sB' % 30 => '30B'
	#struct.pack('30B', *tmp ) has 30 bytes pack binary data.
	#'%dB' % w*h has the same meaning with %d.
	data = struct.pack('%sB' % w*h, *tmp)	#writing newly cropped w*h bytes from *tmp to rid file
	f.write(data)

	#
	f.close()

#write row,col,size(diamter) to list.txt
def export(im, r, c, s, folder, id, list):
	#im.shape=row,col
	nrows = im.shape[0]
	ncols = im.shape[1]
	print('+nrows=', nrows, '+ncols=', ncols)
	# the cropped area is bounded by +-0.75s around center (r,c) which is bigger than
	# face bounding box, +-0.5S, around center <r,c>
	r0 = max(int(r - 0.75*s), 0); r1 = min(r + 0.75*s, nrows)
	c0 = max(int(c - 0.75*s), 0); c1 = min(c + 0.75*s, ncols)
	#the new cropped image containing the original face bounding box.
	im = im[r0:r1, c0:c1]
	#the new dim of the cropped image
	nrows = im.shape[0]
	ncols = im.shape[1]
	print('>nrows=', nrows, '>ncols=', ncols)
	#the new bounding box center in the new cropped image
	r = r - r0
	c = c - c0
	print('>r=',r,' c=',c)
	# resize, if needed
	maxwsize = 192.0
	wsize = max(nrows, ncols)

	ratio = maxwsize/wsize

	if ratio<1.0:#resize the pic because it's > maxwsize
		im = numpy.asarray( Image.fromarray(im).resize((int(ratio*ncols), int(ratio*nrows))) )

		r = ratio*r
		c = ratio*c
		s = ratio*s
	print('resize im= ', im.shape)
	#write "face[id].rid" to dst list.txt
	list.write(id + '.rid\n')

	#creating 7 randomized face bounding boxes from the original bounding box
	nrands = 7;
	#write 7 to dst list.txt
	list.write('\t%d\n' % nrands)

	for i in range(0, nrands):
		#uniformly randomize size (diameter of the bounding box) ratio 0.9-1.1, 10% random
		stmp = s*random.uniform(0.9, 1.1)
		#uniformly randomize row and column, ratio +-5% random
		rtmp = r + s*random.uniform(-0.05, 0.05)
		ctmp = c + s*random.uniform(-0.05, 0.05)

		#plot the 7 bouding boxs
		if plot:
			matplotlib.pyplot.cla()

			matplotlib.pyplot.plot([ctmp-stmp/2, ctmp+stmp/2], [rtmp-stmp/2, rtmp-stmp/2], 'b', linewidth=3)
			matplotlib.pyplot.plot([ctmp+stmp/2, ctmp+stmp/2], [rtmp-stmp/2, rtmp+stmp/2], 'b', linewidth=3)
			matplotlib.pyplot.plot([ctmp+stmp/2, ctmp-stmp/2], [rtmp+stmp/2, rtmp+stmp/2], 'b', linewidth=3)
			matplotlib.pyplot.plot([ctmp-stmp/2, ctmp-stmp/2], [rtmp+stmp/2, rtmp-stmp/2], 'b', linewidth=3)

			matplotlib.pyplot.imshow(im, cmap=matplotlib.cm.Greys_r)

			matplotlib.pyplot.draw()

#			response = input('Press Enter to continue...')
		#write the randomized row, column and size (diameter of bounding box)
		list.write('\t%f %f %f\n' % (rtmp, ctmp, stmp))

	list.write('\n')
	list.flush()

	#write im to dstfolder/facennn.id
	saveasrid(im, folder + '/' + id + '.rid')

def exportmirrored(im, r, c, s, folder, id, list):
	#
	# exploit mirror symmetry of the face
	#

	# flip image
	im = numpy.asarray(ImageOps.mirror(Image.fromarray(im)))

	# flip column coordinate of the object
	c = im.shape[1] - c

	# export
	export(im, r, c, s, folder, id, list)

# image list:
#file0000000000003987.jpg
#file000000000000nnnn.jpg.....
imlist = open(srcfolder + '/Subsets/GENKI-SZSL/GENKI-SZSL_Images.txt', 'r').readlines()

#Each image has three labels: X Position, Y Position, and Size.
#The position labels specify the (x,y) coordinates of the center of a square box surrounding the face.
#The Size label specifies the diameter of the face box.
#The translated object sample is specified by three coordinates (row, column and size; all in pixels)
#rs, an array of: row# of a face bounding box center in the labels
#cs, an array of : col# of a face bounding box center in the labels
#ss, an array of : the diameter of the face box in the labels
rs = [float(line.split()[1]) for line in open(srcfolder+'/Subsets/GENKI-SZSL/GENKI-SZSL_labels.txt', 'r').readlines()]
cs = [float(line.split()[0]) for line in open(srcfolder+'/Subsets/GENKI-SZSL/GENKI-SZSL_labels.txt', 'r').readlines()]
ss = [float(line.split()[2]) for line in open(srcfolder+'/Subsets/GENKI-SZSL/GENKI-SZSL_labels.txt', 'r').readlines()]

#
list = open(dstfolder + '/list.txt', 'w')

n = 0
print('file#=',len(rs))
for i in range(0, len(rs)):
	# image path to /home/thomas/cv/dataset/GENKI/GENKI-R2009a/files/file000000000000nnnn.jpg
	path = srcfolder + '/files/' + imlist[i].strip()
	print(i,':path=',path)
	#get its the bouding facebox, row, column, size(diameter of a bounding box)
	r = rs[i]
	c = cs[i]
	s = ss[i]
	print('r=',r, ' c=',c, 's= ',s);
	#
	try:
		#print('try1')
		#color to grey scale
		#imtmp=Image.open(path)
		#open the jpeg file and covert it from RGB to 8 bit grey scale by convert('L')
		im = Image.open(path).convert('L')
		#im = imtmp.convert('L')
		#print('try2',im)
	except:
		continue
		#print('conti')

	#image im is "width*height" to array im : "row * cols"
	im = numpy.asarray(im)

	#file name as faceid=facennnn.rid in the list.txt
	#write facennnn.rid and parameters to list.txt
	id = 'face' + str(n)
	print('id =', id)
	export(im, r, c, s, dstfolder, id, list)
	n = n+1

	# faces are symmetric and we exploit this here
	#write facennnn+1.rid and parameters to list.txt
	#list.txt has double entries/lines of the GENKI-SZSL_Images.txt,
	#because it has the original and the mirror.
	id = 'face' + str(n)
	exportmirrored(im, r, c, s, dstfolder, id, list)
	n = n+1

#
list.close()
