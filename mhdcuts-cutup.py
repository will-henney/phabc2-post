from PIL import Image
import sys, glob

try:
    runid = sys.argv[1]
except:
    print "Usage: %s RUNID" % sys.argv[0]

filepattern = "hsv-xtd-bbb-cuts-%s-*0.png" % (runid)
filelist = glob.glob(filepattern)

varids = ["bbb", "xtd"]
cutids = ["xy", "xz", "yz"]

nx, ny = 256, 256
mx, my = 8, 8
for imfile in filelist:
    timeid = imfile[-8:-4]
    print imfile
    image = Image.open(imfile)
    oy = my
    for varid in varids: 
        ox = mx
        for cutid in cutids:
            cutout = image.crop((ox, oy, ox+nx, oy+ny))
            cutout.save("hsv-%s-cut-%s-%s-%s.png" % (varid, cutid, runid, timeid))
            ox += nx + mx
        oy += ny + my 
