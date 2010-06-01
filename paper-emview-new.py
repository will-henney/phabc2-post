"""
Figures that show 3D views of the emissivity cubes

The images should already have been created with makemovie

"""
from PIL import Image
import pyx
import axis3d

pyx.text.set(mode="latex")
pyx.text.preamble(r"\AtBeginDocument{\sffamily}")

viewlist = (
#     (310, 310),
    (350, 350),
    (50, 50),
    (10, 80),
    )

# tlist = [20, 40, 50, 60, 70, 80]
tlist = [20, 50, 80]
tlist.reverse()

modelid = "s80-255"
movieid = "evo"
paperid = modelid[:3].upper() + 'H'

# size of images on page
imwidth = 10.0
imheight = imwidth/2
margin = 0.0
x, y = 0, 0

c = pyx.canvas.canvas()
for theta, phi in viewlist:
    for time in tlist:
	# Each view is a row
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

c.writePDFfile("paper-emview-new")
	
	
