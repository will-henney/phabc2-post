from myfitsutils import *
from PIL import Image, ImageDraw
import sys
import colorsys
import numpy as N

def hsv_xtd_image(prefix, it, ix, axis='x'):
    imd = FitsImage(datadir+prefix+'-dd%4.4i.fits' % it, 
		    takelog=1, fmin=5.e-23, fmax=4.e-19, gamma=1.0, 
		    icut=ix, cutaxis=axis)
    imt = FitsImage(datadir+prefix+'-te%4.4i.fits' % it, 
		    takelog=0, fmin=0, fmax=1.e4, icut=ix, cutaxis=axis)
    imx = FitsImage(datadir+prefix+'-xi%4.4i.fits' % it, 
		    takelog=0, fmax=-0.4, fmin=1.05, icut=ix, cutaxis=axis)
    
    # hue is ion frac
    hue = imx.getdata()
    # saturation is temperature
    sat = imt.getdata()
    # value is density
    val = imd.getdata()

    print "XTD %s %s: H = [%3i, %3i] S = [%3i, %3i] V = [%3i, %3i] " \
	% (prefix, axis, N.min(hue), N.max(hue), N.min(sat), N.max(sat), 
	   N.min(val), N.max(val))

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
    return imrgb.transpose(Image.FLIP_TOP_BOTTOM)



def hsv_bbb_image(prefix, it, ix, axis='x'):
    # ax1, ax2 are in the plane of the image
    # ax3 is the perpendicular direction
    if axis=='z':
	ax1 = 'x'
	ax2 = 'y'
	ax3 = 'z'
    elif axis=='y':
	ax1 = 'x'
	ax2 = 'z'
	ax3 = 'y'
    elif axis=='x':
	ax1 = 'y'
	ax2 = 'z'
	ax3 = 'x'
    b1 = FitsData(datadir+prefix+'-b%s%4.4i.fits' % (ax1, it), 
		  icut=ix, cutaxis=axis).fitsdata
    b2 = FitsData(datadir+prefix+'-b%s%4.4i.fits' % (ax2, it), 
		  icut=ix, cutaxis=axis).fitsdata
    b3 = FitsData(datadir+prefix+'-b%s%4.4i.fits' % (ax3, it), 
		  icut=ix, cutaxis=axis).fitsdata
    # magnitude of field
    bb = N.sqrt(b1**2 + b2**2 + b3**2)
    # field angle in plane
    phase = 0.0 #1.0*N.pi
    boost = 2.0			# exaggerate the field angle
    phi = N.remainder(phase + boost*N.arctan2(b1, b2), 
		      2.0*N.pi)/(2.0*N.pi) # should be in range [0, 1]
    # field angle out of plane
    dip = 1.0 - N.abs(b3)/bb # should be in range [0, 1]
    
    b0 = 2.e-4 		# normalization of field
    bgamma = 2.0
    
    # hue is in-plane angle
    hue = phi.flatten()
    # saturation is out-of-plane angle
    sat = dip.flatten()
    # value is field strength
    val = (bb.flatten()/b0)**(1.0/bgamma)

    print "BB %s %s: H = [%.3f, %.3f] S = [%.3f, %.3f] V = [%.3f, %.3f] " \
	% (prefix, axis, hue.min(), hue.max(), sat.min(), sat.max(), 
	   val.min(), val.max())

    # list to contain channel values
    rgb = []
    # convert hsv -> rgb
    for h, s, v in zip(hue.tolist(), sat.tolist(), val.tolist()):
	r, g, b = colorsys.hsv_to_rgb(h, s, v)
	rgb.append( (int(255*r), int(255*g), int(255*b)) )

    # new image for color composite
    imrgb = Image.new('RGB', (bb.shape[1], bb.shape[0]))
    imrgb.putdata(rgb)
    return imrgb.transpose(Image.FLIP_TOP_BOTTOM)

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
    # find the value of y at the peak of globule
    dd = pyfits.open(datadir+runid+'-dd%4.4i.fits' % (i))[0].data
    xn = pyfits.open(datadir+runid+'-xn%4.4i.fits' % (i))[0].data
    dn = (dd*xn).sum(axis=0).sum(axis=-1) # sum over x & z
    jypeak = dn.argmax()
    print 'Y peak is ', jypeak + 1

    # make the four images
    imxy = hsv_xtd_image(runid, i, zcut, axis='z')
    nx, ny = imxy.size
    imxz = hsv_xtd_image(runid, i, jypeak, axis='y')
    imbxy = hsv_bbb_image(runid, i, zcut, axis='z')
    imbxz = hsv_bbb_image(runid, i, jypeak, axis='y')

    ex = 10 			# letter height

    # put the four cuts together
    im4 = Image.new(imxy.mode, (2*nx, 2*ny + 5*ex), color="black")
    im4.paste(imxy, (0, ny))
    im4.paste(imxz, (nx, ny))
    im4.paste(imbxy, (0, 0))
    im4.paste(imbxz, (nx, 0))

    # add labels to the image
    d = ImageDraw.Draw(im4)
    d.text((0.5*ex, ny-1.5*ex), 'XY plane (Z = %i)' % zcut, fill="yellow") 
    d.text((nx+0.5*ex, ny-1.5*ex), 'XZ plane (Y = %i)' % (jypeak+1), fill="yellow") 
    d.text((0.5*ex, ny+0.5*ex), 'XY plane (Z = %i)' % zcut, fill="yellow") 
    d.text((nx+0.5*ex, ny+0.5*ex), 'XZ plane (Y = %i)' % (jypeak+1), fill="yellow") 
    d.text((0.5*ex, 0.5*ex), 'B', fill="yellow") 
    d.text((nx+0.5*ex, 0.5*ex), 'B', fill="yellow") 
    d.text((0.5*ex, 2*ny-1.5*ex), 'x, T, rho', fill="yellow") 
    d.text((nx+0.5*ex, 2*ny-1.5*ex), 'x, T, rho', fill="yellow") 

    d.text((0.5*ex, 2*ny+0.5*ex), 
	   'R = 0.2 pc, D = 0.5 pc, t = %3.3i,000 years' % (i), fill="white") 
    d.text((0.5*ex, 2*ny+2*ex), 'n = 1.0e4 pcc, B = (0.13, 0.13, 0.00) mG', 
	   fill="white") 
    d.text((0.5*ex, 2*ny+3.5*ex), 'Henney, Arthur, de Colle, Mellema (2008)', 
	   fill="white") 
    im4.save('hsv-xtd-bbb-cuts-%s-%4.4i.png' % (runid, i))
    

