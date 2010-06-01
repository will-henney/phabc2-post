"""
Figures that show 3D views of the emissivity cubes

The images should already have been created with makemovie

"""
from PIL import Image
import pyx
import axis3d

pyx.text.set(mode="latex")
pyx.text.preamble(r"\AtBeginDocument{\sffamily}")

theta, phi = (350, 350)

timelist = [30, 60, 90, 120]
timelist.reverse()
modelidlist = ["s00-127m", "s45-127m"]

movieid = "evo"


# size of images on page
imwidth = 10.0
imheight = imwidth/2
margin = 0.0
x, y = 0, 0

c = pyx.canvas.canvas()
for modelid in modelidlist:
    paperid = modelid[:3].upper() + 'L'
    for time in timelist:
	angid = '%+3.3i%+3.3i' % (theta, phi)
	pngfile = "%(modelid)s/rgb-NHO-%(modelid)s-%(movieid)s%(angid)s-"\
	    "%(time)4.4i.png" % locals()
	im = Image.open(pngfile)
	bm = pyx.bitmap.bitmap(x, y, im, width=imwidth)
	c.insert(bm)
	ax = axis3d.axis3d(theta, phi, distance=200, figsize=0.5)
	c.insert(ax, [pyx.trafo.translate(x + (imwidth - 2.0) - ax.ox, y + 1.0 - ax.oy)])
	c.text(x + 0.2, y + 0.2, 
	       r"\boldmath\bfseries %(paperid)s, t = %(time)i,000 years, \(\theta = %(theta)i^\circ\), \(\phi = %(phi)i^\circ\)" % locals(), 
	       [pyx.color.rgb.white])
	y += imheight + margin
    y = 0.0
    x += imwidth + margin

c.writePDFfile("paper-emview-S00L-S45L-new")
	
	
