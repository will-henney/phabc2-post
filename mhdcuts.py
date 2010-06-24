from myfitsutils import *
from PIL import Image
import sys
import colorsys
import numpy as N

def colorimage(prefix, it, ix, axis='x'):
    imd = FitsImage(datadir+prefix+'-dd%4.4i.fits' % it, 
		    takelog=1, fmin=1.e-22, fmax=2.e-19, gamma=1.0, 
		    icut=ix, cutaxis=axis)
    imd.transpose(Image.FLIP_TOP_BOTTOM)
    imt = FitsImage(datadir+prefix+'-te%4.4i.fits' % it, 
		    takelog=0, fmin=0, fmax=1.e4, icut=ix, cutaxis=axis)
    imt.transpose(Image.FLIP_TOP_BOTTOM)
    imx = FitsImage(datadir+prefix+'-xn%4.4i.fits' % it, 
		    takelog=0, fmin=-0.05, fmax=1.4, icut=ix, cutaxis=axis)
    imx.transpose(Image.FLIP_TOP_BOTTOM)
    
#     print 'Read FITS files'
    
    # hue is ion frac
    hue = imx.getdata()
    # saturation is temperature
    sat = imt.getdata()
    # value is density
    val = imd.getdata()

    print "Hue range: ", N.min(hue), N.max(hue)
    print "Sat range: ", N.min(sat), N.max(sat)
    print "Val range: ", N.min(val), N.max(val)

    # list to contain channel values
    rgb = []
    # convert hsv -> rgb
    for h, s, v in zip(hue, sat, val):
	r, g, b = colorsys.hsv_to_rgb(float(h)/255, 
				      float(s)/255, 
				      float(v)/255)
	rgb.append( (int(255*r), int(255*g), int(255*b)) )

#     print 'Converted to HSV->RGB'

    # new image for color composite
    imrgb = Image.new('RGB', imd.size)
    imrgb.putdata(rgb)
    return imrgb


datadir = './'

if len(sys.argv) != 6:
    print "Usage: %s RUNID ZCUT TMIN TMAX TSTEP" % sys.argv[0]
    exit

runid = sys.argv[1]
zcut = int(sys.argv[2])
tmin = int(sys.argv[3])
tmax = int(sys.argv[4])
tstep = int(sys.argv[5])

for i in range(tmin, tmax+1, tstep):
    try:
	im = colorimage(runid, i, zcut, axis='z')
	im.save('mhdcuts-%s-%4.4i.png' % (runid, i))
    except:
	pass


