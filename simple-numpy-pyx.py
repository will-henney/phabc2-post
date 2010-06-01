"""
First stab at redoing the appendix plots using PyX

09 Apr 2009
"""

from heatcool import heat, cool
import numpy as N
import pyx

a = N.linspace(0.0, 1.0)

pyx.text.set(mode="latex")

g = pyx.graph.graphxy(width=10, 
		      x=pyx.graph.axis.lin(title=r"\(a\)"), 
		      y=pyx.graph.axis.lin(title=r"\(f(a)\)"),
		      key=pyx.graph.key.key(pos="tl"))

lines = []
for f, s in [
    [3 * a**2 * (1.0-a), r"\(3 a^2 (1 - a)\)"],
    [0.5*N.sin(N.pi*a**3), r"\(\frac12 \sin \pi a^3\)"],
    [N.log(1.+a), r"\(\ln (1 + a)\)"],
    ]: 
    lines.append(pyx.graph.data.values(x=a, y=f, title=s))

g.plot( lines, [pyx.graph.style.line()] )

basename = __file__.split('.')[0]

g.writePDFfile("%s" % (basename))


