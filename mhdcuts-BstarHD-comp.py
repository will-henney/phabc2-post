from PIL import Image
ntimes = 135
itimes = range(1, ntimes+1)
megayears = [float(i)/100.0 for i in itimes]

filemetapattern = "hsv-xtd-cut-%%s-Bstar-%%s-%4.4i.png"
# Example: hsv-xtd-cut-yz-Bstar-ep-1180.png
filepatterns = [filemetapattern % (10*i) for i in itimes]

runids = ["ep", "HDep"]
cutids = ["xy", "xz", "yz"]

nx, ny = 512, 512                       # Desired size of each pane. Original is scaled up to this
mx, my = 16, 16                         # Margin sizes



##
## Machinery for converting text -> image
##
import pyx, subprocess, hashlib
pyx.text.set(mode="latex")
preambletext = r"""
\usepackage{color}
\usepackage[varg]{txfonts}
%\AtBeginDocument{\pagecolor{black}}
"""
pyx.text.preamble(preambletext)
gscommand = """
gs -sDEVICE=pngalpha -sOutputFile=%s.png -r300 -dBATCH -dNOPAUSE -dSAFER -q -dNOCACHE %s.pdf
"""
tmpfilepattern = "/tmp/mhdcuts-stitchup-%s"

def text2image(text):
    """
    Converts argument (interpreted as LaTeX) into an antialiased PNG image with alpha channel

    Implements simple caching of images to save on LaTeX and gs processing
    """
    # Temporary filename is unique for given combo of text and other relevant stuff
    tmpprefix = tmpfilepattern % (hashlib.md5(text+preambletext+gscommand+textwrapper).hexdigest())
    try:
        # use pre-generated temporary PNG file if possible
        im = Image.open("%s.png" % (tmpprefix))
    except IOError:
        c = pyx.canvas.canvas()
        c.text(0, 0, textwrapper % (text), [pyx.trafo.scale(1.0)])
        c.writePDFfile(tmpprefix)
        subprocess.Popen(gscommand % (tmpprefix, tmpprefix), shell=True).wait()
        im = Image.open("%s.png" % (tmpprefix))
    return im

runtexts = [r"MHD simulation", r"Pure HD simulation"] 
cuttexts = [r"\(xy\) plane", r"\(xz\) plane", r"\(yz\) plane"]
textwrapper = r"{\color{white}\raggedright\bfseries\boldmath %s}"
runtextimages = [text2image(runtext) for runtext in runtexts]
cuttextimages = [text2image(cuttext) for cuttext in cuttexts]
cuttext_height = max([im.size[1] for im in cuttextimages])
##
## Stitch panels together
##
keyimages = [Image.open("mhdcuts-Bkey-%2.2i.png" % (i)) for i in [10, 5]]

for filepattern, itime, megayear in zip(filepatterns, itimes, megayears):
    print filepattern
    rows = [[Image.open(filepattern % (cutid, runid)).resize((nx,ny)) 
             for cutid in cutids] for runid in runids]
    image = Image.new('RGB', (3*nx + 4*mx, 
                              2*ny + 3*my + 2*cuttext_height + mx))
    oy = my
    ox = mx
    for j, row in enumerate(rows):
        ox = mx
        for pane in row:
            image.paste(pane, (ox, oy))
            # image.paste(runtextimages[j], (ox + mx, oy + my), runtextimages[j])
            ox += nx + mx
        oy += ny + my
    ox = mx
    for cuttextimage in cuttextimages:
        dx = (nx - cuttextimage.size[0])/2
        image.paste(cuttextimage, (ox + dx, oy), cuttextimage)
        ox += nx + mx
    oy += cuttext_height + mx
    timetext = r"\(t = %.2f\) Myr" % (megayear)
    timetextimage = text2image(timetext)
    ox = mx
    image.paste(timetextimage, (ox, oy), timetextimage)

    # oy -= 2*mx
    # ox = mx + mx/2 + nx - keyimages[0].size[0]/2
    # image.paste(keyimages[0], (ox, oy), keyimages[0])
    # ox += nx + mx
    # image.paste(keyimages[1], (ox, oy), keyimages[1])

    image.save("mhdcuts-BstarHD-comp-%3.3i.png" % (itime))

