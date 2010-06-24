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
    itime = int(sys.argv[2])
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME [THETA PHI BRIGHTMAX]" % execname
    exit

try: 
    th0 = float(sys.argv[3])
    ph0 = float(sys.argv[4])
except IndexError, ValueError:
    th0, ph0 = 0.0, 0.0
# Optional brightness arguments
try: 
    brightmax0 = float(sys.argv[5])
except IndexError, ValueError:
    print "No brightmax given - using default value"
    brightmax0 = 1.e6

Movie.emtypes = ["S26731", "N26584", "Halpha"]
Movie.bandscales = [0.15, 0.3, 0.8]
movieid="tumble-%4.4i" % itime
movie = Movie(runid=runid, movieid=movieid, 
	      time=itime, nframes=72, verbose=1, force=0)
movie.imageprefix = "rgb-SNH-%s-%s" % (runid, movieid)
movie.brightmaxfunc = brightmax
movie.camera.set_angles(th0, ph0)
movie.boxsize = 4.0
movie.camera.set_steps(5.0, 5.0)
movie.makemovie()
