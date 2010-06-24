"""
Make 2x2 mosaic of visualization images for Bstar:
         +---------------+---------------+
         |               |               |
         |   Bstar-ep    |  Bstar-HDep   |
         |               |               |
         |     SNH       |     SNH       |
         |               |               |
         +---------------+---------------+
         |               |               |
         |   Bstar-ep    |  Bstar-HDep   |
         |               |               |
         |     CPF       |     CPF       |
         |               |               |
         +---------------+---------------+
"""

from PIL import Image

ntimes = 135
itimes = range(1, ntimes+1)
ntumbles = 72
itumbles = range(1, ntumbles+1)

isequence = itimes + itumbles
movietypes = ["evo"]*ntimes + ["tumble"]*ntumbles
thetas = [350]*ntimes + [(350 + (i-1)*5) % 360 for i in itumbles] 
phis = thetas
movienums = ["+350+350"]*ntimes + ["-1350"]*ntumbles
megayears = [float(i)/100.0 for i in itimes + [ntimes]*ntumbles]
filemetapattern = "rgb-%%s-Bstar-%%s-%s%s-%4.4i.png"
# Example: rgb-CPF-Bstar-ep-evo+350+350-0128.png
filepatterns = [filemetapattern % (mtype, mnum, i) 
                for mtype, mnum, i in zip(movietypes, movienums, isequence)]  

chanids =  ["SNH", "CPF"]
runids = ["ep", "HDep"]
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
tmpfilepattern = "/tmp/stitchup-%s"

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
        c.text(0, 0, textwrapper % (text))
        c.writePDFfile(tmpprefix)
        subprocess.Popen(gscommand % (tmpprefix, tmpprefix), shell=True).wait()
        im = Image.open("%s.png" % (tmpprefix))
    return im

##
## Make text for each channel 
##
chantexts = [
    r"""
Hubble WFPC2/WFC3:\\
\textcolor{red}{\(\bullet\) [S\,\textsc{ii}] 6731~\AA}\\
\textcolor{green}{\(\bullet\) [N\,\textsc{ii}] 6584~\AA}\\
\textcolor{blue}{\(\bullet\) H\(\alpha\) 6563\AA}\\
    """,
    r"""
Herschel/Spitzer/VLA:\\
\textcolor{red}{\(\bullet\) Far-IR Dust}\\
\textcolor{green}{\(\bullet\) Mid-IR PAH}\\
\textcolor{blue}{\(\bullet\) 6\,cm Free-Free}\\
    """,
    ]
textwrapper = r"\parbox{4cm}{\color{white}\raggedright\bfseries\boldmath %s}"
chantextimages = [text2image(chantext) for chantext in chantexts]
chantext_width = max([im.size[0] for im in chantextimages])


## 
## Text for HD/MHD
## 
runtexts = [r"MHD simulation", r"Pure HD simulation"]
textwrapper = r"\framebox{\color{white}\raggedright\bfseries\boldmath %s}"
runtextimages = [text2image(runtext) for runtext in runtexts]
runtext_height = max([im.size[1] for im in runtextimages])

##
## Stitch panels together
##
jcount = 1
for filepattern, theta, phi, megayear in zip(filepatterns, thetas, phis, megayears):
    print filepattern
    rows = [[Image.open(filepattern % (chanid, runid)).resize((nx,ny)) 
              for runid in runids] for chanid in chanids]
    image = Image.new('RGB', (2*nx + 3*mx + chantext_width + mx, 
                              2*ny + 3*my + 2*runtext_height + 2*my))
    oy = my
    ox = mx
    for runtextimage in runtextimages:
        image.paste(runtextimage, (ox, oy), runtextimage)
        ox += nx + mx
    oy += runtext_height + my
    for row, chantextim in zip(rows, chantextimages):
        ox = mx
        for pane in row:
            image.paste(pane, (ox, oy))
            ox += nx + mx
        dy = (ny - chantextim.size[1])/2 # center text vertically
        image.paste(chantextim, (ox, oy+dy), chantextim) # third argument is mask to use alpha channel
        oy += ny + my
    #
    # Text for time/angle info
    #
    timetext = r"\(t = %.2f\) Myr \quad \(\theta = %3.3i^\circ\)\quad \(\phi = %3.3i^\circ\)" % (megayear, theta, phi)
    timetextimage = text2image(timetext)
    ox = mx
    image.paste(timetextimage, (ox, oy), timetextimage)

    image.save("movie-Bstar-stitchup-%3.3i.png" % (jcount))
    jcount += 1

