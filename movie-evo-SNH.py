import sys, os
from movierot import Movie

def brightmax(i): 
    """
    Smoothly varying brightness with time so that the movie looks good
    """
    return brightmax0/(1.0+(float(i)/500)**2)


# Parse command line arguments
execname = os.path.split(sys.argv[0])[-1]
try: 
    runid = sys.argv[1]
    it1 = int(sys.argv[2])
    it2 = int(sys.argv[3])
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME1 ITIME2 [IDELTA THETA PHI BRIGHTMAX]" % execname
    exit

# Optional orientation arguments
try: 
    idelta = int(sys.argv[4])
except IndexError, ValueError:
    idelta = 1

try: 
    th0 = float(sys.argv[5])
    ph0 = float(sys.argv[6])
except IndexError, ValueError:
    th0, ph0 = 0.0, 0.0
# Optional brightness arguments
try: 
    brightmax0 = float(sys.argv[7])
except IndexError, ValueError:
    print "No brightmax given - using default value"
    brightmax0 = 1.e6


Movie.emtypes = ["S26731", "N26584", "Halpha"]
Movie.bandscales = [0.15, 0.3, 0.8]
movieid = "evo%+3.3i%+3.3i" % (th0, ph0)
movie = Movie(runid=runid, movieid=movieid, 
	      time=it1, dtime=idelta, frame=it1/idelta, 
              nframes=1+(it2-it1)/idelta, verbose=1, force=0)
movie.imageprefix = "rgb-SNH-%s-%s" % (runid, movieid)
movie.brightmaxfunc = brightmax
movie.camera.set_steps(0.0, 0.0)
movie.camera.set_angles(th0, ph0)
movie.boxsize = 4.0

movie.makemovie()
