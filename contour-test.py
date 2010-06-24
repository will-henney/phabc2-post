"""
Stress test for pyxgraph's pyxplotarray and pyxplotcontour

Problems found with pyxgraph 0.2.5: 

1. pyxplotarray shifts the image 0.5 pixel to lower right
2. pyxplotcontour misses out the topmost row and rightmost column

Will Henney - 31 Dec 2007
"""
import numpy, pyxgraph

nx, ny = 12, 15
# coordinates of pixel centers
x = 0.5 + numpy.arange(nx)[numpy.newaxis,:]
y = 0.5 + numpy.arange(ny)[:,numpy.newaxis]

# pyxplotcontour only works with 1D arrays
xx = 0.5 + numpy.arange(nx)
yy = 0.5 + numpy.arange(ny)

x0, y0 = (float(nx)/2, float(ny)/2) # origin for function
r = numpy.sqrt((x - x0)**2 + (y - y0)**2) # radius from origin

# # Gaussian function
# sigma = float(n)/20		  # width
# z = numpy.exp(-r**2/sigma**2)

# Radius function
z = r

pixmargin = 0

figwidth = 10
figheight = (ny+2*pixmargin)*figwidth/(nx + 2*pixmargin)
g = pyxgraph.pyxgraph(xlimits=(-pixmargin, nx+pixmargin), 
		      ylimits=(-pixmargin, ny+pixmargin),
		      width=figwidth, height=figheight, key=None)
graymap = pyxgraph.ColMapper.ColorMapper("pm3d", pm3d=[3,3,3])
g.pyxplotarray(z, colmap=graymap, compressmode=None,
	       xpos=0.0, ypos=0.0, width=nx, height=ny, graphcoords=True)
g.pyxplotcontour(z, xx, yy, levels=4, colors='color', color='red')
# plot a big green cross where the peak should be
g.pyxplot(data=((x0,0,0,nx,nx),(y0,0,ny,0,ny)), style="points", color="green", ps=5, pt=0)
g.writePDFfile("contour-test")
