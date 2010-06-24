from myfitsutils import *
from PIL import Image
import sys
import pyx
import numpy as N

datadir = './'

if len(sys.argv) != 5:
    print "Usage: %s RUNID TIME AXIS CUT" % sys.argv[0]
    exit

runid = sys.argv[1]
itime = int(sys.argv[2])
axis = sys.argv[3]
icut = int(sys.argv[4])


imd = FitsImage(datadir+prefix+'-dd%4.4i.fits' % it, 
		takelog=1, fmin=2.e-21, fmax=5.e-24, gamma=1.0, 
		icut=ix, cutaxis=axis)
imd.transpose(Image.FLIP_TOP_BOTTOM)
imt = FitsImage(datadir+prefix+'-te%4.4i.fits' % it, 
		takelog=0, fmin=0, fmax=1.e4, icut=ix, cutaxis=axis)
imt.transpose(Image.FLIP_TOP_BOTTOM)
imx = FitsImage(datadir+prefix+'-xn%4.4i.fits' % it, 
		takelog=0, fmin=-0.05, fmax=1.4, icut=ix, cutaxis=axis)
imx.transpose(Image.FLIP_TOP_BOTTOM)

