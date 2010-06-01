"""
Figures that show 3D views of the emissivity cubes

The images should already have been created with makemovie

"""
from PIL import Image
import pyx
import axis3d

pyx.text.set(mode="latex")
pyx.text.preamble(r"\usepackage{mathptmx}\AtBeginDocument{\sffamily}")

viewlist = (
    (5, 5),
    )

tlist = (60, 90)

modelid = "glob13-128x"
movieid = "tumble"


# size of images on page
imwidth = 10.0
imheight = imwidth/2
margin = 0.5
x, y = 0, 0

c = pyx.canvas.canvas()
for time in tlist:
    # Each time is a row
    for theta, phi in viewlist:
	# Each view is a column
	angid = theta/5.0
	pngfile = "rgb-NHO-%(modelid)s-%(movieid)s-"\
	    "%(time)4.4i-%(angid)4.4i.png" % locals()
	im = Image.open(pngfile)
	bm = pyx.bitmap.bitmap(x, y, im, width=imwidth)
	c.insert(bm)
	ax = axis3d.axis3d(theta, phi, distance=200, figsize=0.5)
	c.insert(ax, [pyx.trafo.translate(x + (imwidth - 2.0) - ax.ox, y + 1.0 - ax.oy)])

	x += imwidth + margin
    x = 0.0
    y -= imheight + margin

c.writePDFfile("paper-emview-Z00L")
	
	
