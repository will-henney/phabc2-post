from myfitsutils import *
from PIL import Image, ImageDraw
import sys
import colorsys
import numpy as N

def hsv_xtd_image(prefix, it, ix, axis='x'):
    imd = FitsImage(datadir+prefix+'-dd%4.4i.fits' % it, 
		    takelog=1, fmin=rhomin, fmax=rhomax, gamma=1.0, 
		    icut=ix, cutaxis=axis)
    imt = FitsImage(datadir+prefix+'-te%4.4i.fits' % it, 
		    takelog=1, fmin=100.0, fmax=1.e4, gamma=2.0, icut=ix, cutaxis=axis)
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
    
    b0 = 6.e-6 		# normalization of field
    bgamma = 1.0
    
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

if "krumx" in runid:
    rhomin, rhomax = 2.e-27, 1.e-20
else:
    rhomin, rhomax = 2.e-25, 4.e-21

for i in range(tmin, tmax+1, tstep):
    print "........................................................................"
    print "time: ", i
    print "........................................................................"

    try: 
	# make the six images
	imxy = hsv_xtd_image(runid, i, zcut, axis='z')
	print "Images are %i x %i" % (imxy.size)
	nx, ny = imxy.size
	imxz = hsv_xtd_image(runid, i, zcut, axis='y')
	imyz = hsv_xtd_image(runid, i, zcut, axis='x')
	imbxy = hsv_bbb_image(runid, i, zcut, axis='z')
	imbxz = hsv_bbb_image(runid, i, zcut, axis='y')
	imbyz = hsv_bbb_image(runid, i, zcut, axis='x')
	print "B images are %i x %i" % (imbxy.size)
    except: 
        print "Could not make the images!"
	continue

    ex = 10 			# letter height
    m = 8 # margin

    # put the four cuts together (margins left and right)
    im6 = Image.new(imxy.mode, (3*nx+4*m, 2*ny + 3*m + 6.5*ex), color="black")
    im6.paste(imxy, (m + 0, m + ny + m))
    im6.paste(imxz, (m + nx + m, m + ny + m))
    im6.paste(imyz, (m + 2*(nx + m), m + ny + m))
    im6.paste(imbxy, (m + 0, m + 0))
    im6.paste(imbxz, (m + nx + m, m + 0))
    im6.paste(imbyz, (m + 2*(nx + m), m + 0))

    # add labels to the image
    d = ImageDraw.Draw(im6)
    d.text((m+0.5*ex, m+ny-0.3*ex), 'XY plane', fill="yellow") 
    d.text((m+nx+m+0.5*ex, m+ny-0.3*ex), 'XZ plane', fill="yellow") 
    d.text((m+2*(nx+m)+0.5*ex, m+ny-0.3*ex), 'YZ plane', fill="yellow") 
    d.text((m+0.5*ex, m+0.2*ex), 'B', fill="yellow") 
    d.text((m+nx+m+0.5*ex, m+0.2*ex), 'B', fill="yellow") 
    d.text((m+2*(nx+m)+0.5*ex, m+0.2*ex), 'B', fill="yellow") 
    d.text((m+0.5*ex, 2*(ny+m)-1.5*ex), 'x, T, rho', fill="yellow") 
    d.text((m+nx+m+0.5*ex, 2*(ny+m)-1.5*ex), 'x, T, rho', fill="yellow") 
    d.text((m+2*(nx+m)+0.5*ex, 2*(ny+m)-1.5*ex), 'x, T, rho', fill="yellow") 

    d.text((m, 2*(ny+m)+0.5*ex), 
           'Run %s, t = %3.3i0 kyr' % (runid, i), fill="white") 
    d.text((m, 2*(ny+m)+2*ex), 'Cuts at X, Y, Z = %i' % zcut, 
	   fill="white") 
    d.text((m, 2*(ny+m)+3.5*ex), 'Arthur, Henney, de Colle,', 
	   fill="white") 
    d.text((m, 2*(ny+m)+5*ex), 'Mellema, Vazquez-Semadeni (2010)', 
	   fill="white") 
    im6.save('hsv-xtd-bbb-cuts-%s-%4.4i.png' % (runid, i))
    

