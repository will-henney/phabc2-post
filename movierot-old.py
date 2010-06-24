"""
Make a movie of a rotating view of a datacube at a fixed time
"""
import subprocess, StringIO, os, sys
import myfitsutils as F
from PIL import Image


class Camera(object):
    """
    A camera with (theta, phi) orientation
    """
    theta = 0.0
    phi = 0.0
    def set_angles(self, theta, phi):
	self.theta = theta
	self.phi = phi
    def set_steps(self, dth, dph):
	self.dth = dth
	self.dph = dph
    def __init__(self, theta=0.0, dth=1.0, phi=0.0, dph=1.0):
	# Initial angles
	self.set_angles(theta, phi)
	# Atomic steps for angles
	self.set_steps(dth, dph)
    def advance(self, incr=1):
	"Advance angles by increment"
	self.theta += incr*self.dth
	self.phi += incr*self.dph
	

class Movie(object):
    runid = None
    time = 100
    emtypes = ["N26584", "Halpha", "O35007"]
    camera = Camera()
    datadir = '.'
    boxsize = 1.0 		# z extent in parsecs
    brightmax = 1.e7		# typical maximum brightness
    # relative scaling of the emission bands
    bandscales = [0.28, 1.0, 0.2]
    gamma = 2.0			# default but can be overwritten by makeRGBmap

    def __init__(self, runid, movieid, 
		 time=100, nframes=100, verbose=0, force=0):
	self.runid = runid
	self.movieid = movieid
	self.time = time
	self.nframes = nframes
	self.verbose = verbose
	self.force = force 	# Whether to regenerate all files
	self.frame = 1		# all movies start with frame 1
	self.execdir = os.path.split(sys.argv[0])[0]
	self.imageprefix = "rgb-NHO-%s-%s" % (self.runid, self.movieid)
	self.imagelist = []
	
    def mapprefix(self, emtype):
	"""
	Returns the prefix string for a given map

	E.g., glob12-128xmap-rot+060-015-Halpha0060
	"""
	return "%smap-rot%+3.3i%+3.3i-%s%4.4i" % (
	    self.runid,
	    self.camera.theta,
	    self.camera.phi,
	    emtype,
	    self.time,
	    )

    def makeRGBimage(self, gamma=None):
	"""
	Make a single RGB emision image

	Returns filename of image
	"""
	imagename = self.imageprefix + "-%4.4i.png" % self.frame
	if self.force > 0 or not os.path.isfile(imagename):
	    if self.verbose > 0: print "Creating image file %s" % imagename
	    rfits, gfits, bfits = [self.makefitsmap(emtype) 
				   for emtype in self.emtypes]
	    rmax, gmax, bmax = [scale*self.brightmax for scale in self.bandscales]
	    if gamma is None: gamma = self.gamma
	    image = F.RGB3Image(rfits, gfits, bfits, 
				0, 0, 0, rmax, gmax, bmax, gamma)
	    image = image.transpose(Image.FLIP_TOP_BOTTOM)
	    image.save(imagename)
	else:
	    if self.verbose > 0: print "Image file %s already exists" % imagename
	return imagename


    def makeimages(self):
	framerange = range(self.frame, self.frame+self.nframes)
	self.imagelist = []
	for self.frame in framerange:
	    self.imagelist.append(self.makeRGBimage())
	    self.camera.advance(1)

    def makemovie(self):
	if len(self.imagelist) != self.nframes:
	    self.makeimages()
	encode_exec = "mencoder"
	encode_args = "-ovc lavc -lavcopts vbitrate=5000:vcodec=wmv2 -mf type=png:fps=15"
	file_args = r"-o %s.avi mf://%s\*.png" % (self.imageprefix, self.imageprefix)
	cmd = "%s %s %s" % (encode_exec, encode_args, file_args)
	subprocess.Popen(cmd, shell=True).wait()
	print "Written %s.avi" % (self.imageprefix)

    def makefitsmap(self, emtype):
	"""
	Make a single FITS emission map

	Returns fits file name
	"""
	fitsfilename = self.mapprefix(emtype) + ".fits"
	if self.force > 1 or \
		not os.path.isfile(os.path.join(self.datadir, fitsfilename)):
	    # try and make the fits file if it is not there, or if force flag is set
	    if self.verbose > 0: print "Creating FITS file %s" % fitsfilename
	    if self.verbose > 1:
		output = None 		# print out all shell output
	    else:
		output = subprocess.PIPE # we never read from this
	    cmd = os.path.join(self.execdir, "makerotmap")
	    p = subprocess.Popen([cmd], stdin=subprocess.PIPE, stdout=output, 
				 cwd=self.datadir)
	    # we should check here that the emissivity cube exists
	    p.stdin.writelines(
		[str(x) + '\n' for x in [
			self.runid,
			self.time,
			emtype,
			self.boxsize,
			self.camera.theta,
			self.camera.phi,
			]]
		)
	    p.stdin.close()
	    p.wait()
	else:
	    if self.verbose > 0: print "FITS file %s already exists" % fitsfilename

	return fitsfilename


if __name__ == '__main__': 
    # Test the machinery
    movie = Movie("glob12-128x", "testmovie100", 
		  time=100, nframes=70, verbose=2, force=0)
    movie.brightmax = 5.e7
    movie.camera.set_steps(5.0, 5.0)
    movie.makemovie()
