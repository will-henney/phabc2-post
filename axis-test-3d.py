import pyx, subprocess

pyx.text.set(mode="latex")
pyx.text.preamble(r"\usepackage{mathptmx}\AtBeginDocument{\sffamily}")
figsize = 0.7
pyx.unit.set(xscale=1)		# scale fonts

def pdf2jpg(filestem):
    "Converts <filestem>.pdf to <filestem>.jpg using Ghostscript"
    gs_res = 300
    gs_opts = "-dBATCH -dNOPAUSE -dSAFER -q -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -dNOCACHE -r%i -sDEVICE=jpeg" % (gs_res)
    gs_cmd = "gs " + gs_opts + " -sOutputFile=%s.jpg %s.pdf"
    subprocess.Popen(gs_cmd % (filestem, filestem), shell=True).wait()

DD = 5

for i in range(360./DD):
    # CAMERA POSITION
    distance = 10
    # Viewing angle in frame of the simulations
    # -ZZ is the LOS vector, (XX, YY) are in POS 
    TH = DD*i			# angle between X and XX
    PH = TH				# angle between Y and YY

    # Transform to PyX viewing angles
    phi = (180 + PH - 90) % 360.0
    theta = (TH) % 360.0

    # GRAPH CANVAS
    blankaxis = pyx.graph.axis.linear(min=0, max=1, painter=None)
    projector = pyx.graph.graphxyz.central(distance, phi, theta)
    g = pyx.graph.graphxyz(size=figsize, xscale=1, yscale=1, zscale=1, projector=projector, x=blankaxis, y=blankaxis, z=blankaxis)

    vg = g.vgeodesic		# short alias to save typing

    # ARROWS
    # arrow shape 
    aa = 1.0			# arrow length
    a = 0.2			# arrow head width
    b = 0.7			# stalk/head boundary
    c = 0.07			# arrow stalk width
    d = 0.5*c			# stalk start

    # z-axis arrow
    p3 = vg(-c, 0, d, c, 0, d) << vg(c, 0, d, c, 0, b) << vg(c, 0, b, a, 0, b) << vg(a, 0, b, 0, 0, aa) << vg(0, 0, aa, -a, 0, b) << vg(-a, 0, b, -c, 0, b)
    p3.append(pyx.path.closepath())

    # Next two axes need flipping
#     aa = -aa
#     b = -b
#     d = -d

    # x-axis arrow
    p1 = vg(d, -c, 0, d, c, 0) << vg(d, c, 0, b, c, 0) << vg(b, c, 0, b, a, 0) << vg(b, a, 0, aa, 0, 0) << vg(aa, 0, 0, b, -a, 0) << vg(b, -a, 0, b, -c, 0)
    p1.append(pyx.path.closepath())

    # y-axis arrow
    p2 = vg(0, d, -c, 0, d, c) << vg(0, d, c, 0, b, c) << vg(0, b, c, 0, b, a) << vg(0, b, a, 0, aa, 0) << vg(0, aa, 0, 0, b, -a) << vg(0, b, -a, 0, b, -c)
    p2.append(pyx.path.closepath())


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
#     r1 = [-rr, 0, 0]
#     r2 = [0, -rr, 0]
    r1 = [rr, 0, 0]
    r2 = [0, rr, 0]
    r3 = [0, 0, rr]
    aligns = [pyx.text.halign.center, pyx.text.valign.middle]

    # DRAW ON CANVAS
    pcrslist = zip([p1, p2, p3], [c1, c2, c3], [cc1, cc2, cc3], [r1, r2, r3], [s1, s2, s3])

    # Make sure we start from the back
    def compare_depth(a, b):
	"Compare LOS depth of r"
	return cmp(-g.vzindex(*a[3]), -g.vzindex(*b[3]))
    pcrslist.sort(compare_depth)

    # We want the X axis (= PyX z axis) to be horizontal, not vertical
    ox, oy = g.pos(0, 0, 0)
    print "Origin at (%.2f, %.2f)" % (ox, oy)
    rotate = pyx.trafo.rotate(-90, x=ox, y=oy)

    for p, c, cc, r, s in pcrslist:
	g.fill(p, [c, alpha, rotate, pyx.deco.stroked([cc])])
	x, y = rotate.apply(*g.vpos(*r))
	# Add black border round axis label to make it stand out
	g.text(x, y, r"\textbf{%s}" % (s), aligns + [pyx.color.rgb.black, pyx.trafo.scale(1.25,1.05)])
	g.text(x, y, r"\textbf{%s}" % (s), aligns + [cc])

    # add to fixed size canvas
    canvas = pyx.canvas.canvas()
    bgw, bgh = 10*figsize, 6*figsize
    bgbox = pyx.path.rect(0, 0, bgw, bgh)
    canvas.fill(bgbox, [pyx.color.rgb.black])
    canvas.insert(g, [pyx.trafo.translate(0.5*bgw-ox, 0.5*bgh-oy)])
    canvas.text(bgw - 6*pyx.unit.x_pt, 6*pyx.unit.x_pt, r'\(\Theta = %3.3i^\circ\), \(\Phi = %3.3i^\circ\)' % (int(TH), int(PH)), [pyx.text.halign.right, pyx.text.valign.bottom, pyx.color.rgb.white])

    outfile = "axis-test-3d-%3.3i" % (int(TH))
    canvas.writetofile(outfile + ".pdf")
    print "written %s.pdf" % outfile
    pdf2jpg(outfile)
    print "written %s.jpg" % outfile
