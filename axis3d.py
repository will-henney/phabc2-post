import pyx, subprocess

class axis3d(pyx.graph.graphxyz):
    
    def __init__(self, Theta, Phi=None, distance=10, figsize=0.7):
	if Phi is None: Phi = Theta
	# Transform to PyX viewing angles
	phi = (180 + Phi - 90) % 360.0
	theta = (Theta) % 360.0

	# GRAPH CANVAS
	blankaxis = pyx.graph.axis.linear(min=0, max=1, painter=None)
	projector = pyx.graph.graphxyz.central(distance, phi, theta)
	pyx.graph.graphxyz.__init__(self, size=figsize, xscale=1, yscale=1, zscale=1, projector=projector, x=blankaxis, y=blankaxis, z=blankaxis)

	vg = self.vgeodesic		# short alias to save typing

	# ARROWS
	# arrow shape 
	aa = 1.0			# arrow length
	a = 0.2			# arrow head width
	b = 0.7			# stalk/head boundary
	c = 0.07			# arrow stalk width
	d = 0.5*c			# stalk start

	# x-axis arrow
	p1 = vg(d, -c, 0, d, c, 0) << vg(d, c, 0, b, c, 0) << vg(b, c, 0, b, a, 0) << vg(b, a, 0, aa, 0, 0) << vg(aa, 0, 0, b, -a, 0) << vg(b, -a, 0, b, -c, 0)
	p1.append(pyx.path.closepath())

	# y-axis arrow
	p2 = vg(0, d, -c, 0, d, c) << vg(0, d, c, 0, b, c) << vg(0, b, c, 0, b, a) << vg(0, b, a, 0, aa, 0) << vg(0, aa, 0, 0, b, -a) << vg(0, b, -a, 0, b, -c)
	p2.append(pyx.path.closepath())

	# z-axis arrow
	p3 = vg(-c, 0, d, c, 0, d) << vg(c, 0, d, c, 0, b) << vg(c, 0, b, a, 0, b) << vg(a, 0, b, 0, 0, aa) << vg(0, 0, aa, -a, 0, b) << vg(-a, 0, b, -c, 0, b)
	p3.append(pyx.path.closepath())

	# (x, y, z) are PyX axes
	# (X, Y, Z) are simulation axes
	# Permutation between the two frames: (X, Y, Z) = (z, x, y)
	# Labels 1, 2, 3 refer to x, y, z


	# COLORS
	alpha = pyx.color.transparency(0.3)
	ca = 1				# intensity of main channel
	cb = 0.5			# intensity of other channels
	c3 = pyx.color.rgb(ca, cb, cb)	# X is red
	c1 = pyx.color.rgb(cb, ca, cb)	# Y is green
	c2 = pyx.color.rgb(cb, cb, ca)	# Z is blue
	cc3 = pyx.color.rgb(ca, 0, 0)
	cc1 = pyx.color.rgb(0, ca, 0)
	cc2 = pyx.color.rgb(0, 0, ca)


	# AXIS LABELS
	s3 = 'X'
	s1 = 'Y'
	s2 = 'Z'
	rr = 1.1 			# position of label center along each axis
	r1 = [rr, 0, 0]
	r2 = [0, rr, 0]
	r3 = [0, 0, rr]
	aligns = [pyx.text.halign.center, pyx.text.valign.middle]

	# DRAW ON CANVAS
	pcrslist = zip([p1, p2, p3], [c1, c2, c3], [cc1, cc2, cc3], [r1, r2, r3], [s1, s2, s3])

	# Make sure we start from the back
	def compare_depth(a, b):
	    "Compare LOS depth of r"
	    return cmp(-self.vzindex(*a[3]), -self.vzindex(*b[3]))
	pcrslist.sort(compare_depth)

	# We want the X axis (= PyX z axis) to be horizontal, not vertical
	self.ox, self.oy = self.pos(0, 0, 0)
	rotate = pyx.trafo.rotate(-90, x=self.ox, y=self.oy)

	shadowcolor = pyx.color.rgb.black
	shadowopacity = pyx.color.transparency(0.5)
	shadowshift = pyx.trafo.translate(1*pyx.unit.x_pt, -1*pyx.unit.x_pt)
	shadowattrs = [shadowcolor, shadowopacity, shadowshift]
	arrowstrokeattrs = [pyx.style.linewidth.Thick, pyx.style.linejoin.round]
	for p, c, cc, r, s in pcrslist:
	    self.fill(p, [c, alpha, rotate, pyx.deco.stroked([cc] + arrowstrokeattrs)])
	    x, y = rotate.apply(*self.vpos(*r))
	    axistext = r"\textbf{%s}" % (s)
	    # Add shadow to text to make it stand out
	    self.text(x, y, axistext, aligns + shadowattrs)
	    self.text(x, y, axistext, aligns + [cc])





#     # add to fixed size canvas
#     canvas = pyx.canvas.canvas()
#     bgw, bgh = 10*figsize, 6*figsize
#     bgbox = pyx.path.rect(0, 0, bgw, bgh)
#     canvas.fill(bgbox, [pyx.color.rgb.black])
#     canvas.insert(g, [pyx.trafo.translate(0.5*bgw-ox, 0.5*bgh-oy)])
#     canvas.text(bgw - 6*pyx.unit.x_pt, 6*pyx.unit.x_pt, r'\(\Theta = %3.3i^\circ\), \(\Phi = %3.3i^\circ\)' % (int(TH), int(PH)), [pyx.text.halign.right, pyx.text.valign.bottom, pyx.color.rgb.white])

#     outfile = "axis-test-3d-%3.3i" % (int(TH))
#     canvas.writetofile(outfile + ".pdf")
#     print "written %s.pdf" % outfile
#     pdf2jpg(outfile)
#     print "written %s.jpg" % outfile
