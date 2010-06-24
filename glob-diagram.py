import pyx
import myfitsutils as F
import numpy as N


runid = 'glob12-256x'
itime = 0
var = 'dd'
filename = '%s-%s%4.4i.fits' % (runid, var, itime)
maxden = 1.0e4*1.3*1.67262158e-24
im = F.FitsImage(filename, 
		 fmin=1.3*maxden, fmax=0.0, 
		 gamma = 0.5,
		 cutaxis='z', icut=128)
imwidth = 10.0
imheight = imwidth/2
gmargin = 0.1

# setup pyx 
pyx.text.set(mode="latex", errordebug=2, texdebug="glob-diagram.log")
pyx.text.preamble(r"""%
\usepackage{color}
\newcommand\B[1]{\colorbox{white}{#1}}
\setlength\fboxsep{1pt}
""")

textattrs = [pyx.text.size(-2), 
#  	     pyx.text.valign.middle, 
	     pyx.text.vshift.mathaxis, 
	     pyx.text.halign.center]

# Make the graph, with margin around simulation box
g = pyx.graph.graphxy(
    width=imwidth*(1+gmargin), height=imheight*(1+2*gmargin),
    x = pyx.graph.axis.linear(title=r"$x$, parsec", 
			      min = 0.0-gmargin, max = 2.0+gmargin),
    y = pyx.graph.axis.linear(title=r"$y$, parsec",
			      min = 0.0-gmargin, max = 1.0+gmargin),
    )
g.finish()

x0, y0 = g.pos(0, 0)		# Origin of simulation box
# Add image of initial density field
g.insert(pyx.bitmap.bitmap(x0, y0, im, width=imwidth))

def cosd(theta):
    return N.cos(N.pi*theta/180)
def sind(theta):
    return N.sin(N.pi*theta/180)

bline = pyx.path.line(x0-imwidth, y0, x0+3*imwidth, y0)
theta = 60
nlines = 10
dy = imheight/nlines/cosd(theta)

# stroke the y = 0.5 axis
g.stroke(g.ygridpath(0.5), [pyx.style.linestyle.dotted])

# The extent of the simulation box
boxpath = pyx.path.rect(x0, y0, imwidth, imheight)
boxclip = pyx.canvas.canvas([pyx.canvas.clip(boxpath)])

# Draw the magnetic field lines, clipped to boxpath
x0, y0 = g.pos(0, 0.5)		# 
for i in range(-2*nlines, 2*nlines):
    x, y = x0, y0+i*dy
    translate = pyx.trafo.translate(x, y)
    rotate = pyx.trafo.rotate(theta, x, y)
    boxclip.stroke(bline, [translate, rotate])

# Insert into graph
g.insert(boxclip)
# And draw the box
g.stroke(boxpath, [pyx.style.linestyle.dashed])

# Label the field lines and mark the angle
x, y = g.pos(0.9/sind(theta), 0.5)
translate = pyx.trafo.translate(x, y)
rotate = pyx.trafo.rotate(theta, 0, 0)
g.text(0, 0, r"\B{B-field}", [rotate, translate] + textattrs)

x, y = x0 + 5*dy/sind(theta), y0
em = 12*pyx.unit.x_pt

arcpath = pyx.path.path(pyx.path.arc(x, y, 0.5, 0, theta))
g.stroke(arcpath, [pyx.deco.earrow.normal])
g.text(x+0.6*em, y+0.3*em, r"$\theta$", textattrs)

# Add the star
x, y = g.pos(0, 0.5)
g.fill(pyx.path.circle(x, y, 0.15))
g.text(x, y - 0.3, r"\B{star}", textattrs)

# Label the globule
x, y = g.pos(0.5, 0.5)
g.text(x, y - 1.5, r"\B{globule}", textattrs)

# Label the globule
x, y = g.pos(1.0, 0.0)
g.text(x, y, r"\B{simulation box}", textattrs)

g.writePDFfile("glob-diagram")
g.writePSfile("glob-diagram")



