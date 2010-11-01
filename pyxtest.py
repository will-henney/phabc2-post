import pyx
pyx.text.set(mode="latex")
pyx.text.preamble("""\usepackage{mathptmx}""")
c = pyx.canvas.canvas()
c.text(0, 0, "Test that all is well with PyX")
c.writePDFfile("pyxtest")
