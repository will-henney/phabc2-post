import sys, os
from movierot import Movie

def brightmax(i): 
    """
    Smoothly varying brightness with time so that the movie looks good
    """
    return brightmax0/(1.0+(float(i)/100)**2)


# Parse command line arguments
execname = os.path.split(sys.argv[0])[-1]
try: 
    runid = sys.argv[1]
    itime = int(sys.argv[2])
except (IndexError, ValueError):
    print "Usage: %s RUNID ITIME [THETA PHI " % execname
    exit

try: 
    th0 = float(sys.argv[3])
    ph0 = float(sys.argv[4])
except (IndexError, ValueError):
    th0, ph0 = 0.0, 0.0

# Optional brightness arguments
# set defaults
brightmax0 = 2.e8
bandscales = Movie.bandscales
try: 
    # 5th arg is scalar: use to set brightmax0
    brightmax0 = float(sys.argv[5])
except IndexError:
    pass
except ValueError:
    try:
        # 5th arg is triplet: use to set bandscales
        bandscales = [float(x) for x in sys.argv[5].split()]
    except ValueError:
        # otherwise, just use defaults
        pass

Movie.bandscales = bandscales
movieid="tumble-%4.4i" % itime
movie = Movie(runid=runid, movieid=movieid, 
	      time=itime, nframes=72, verbose=1, force=0)
movie.imageprefix = "rgb-NHO-%s-%s" % (runid, movieid)
movie.brightmaxfunc = brightmax
movie.boxsize = 4.0
movie.camera.set_steps(5.0, 5.0)
movie.makemovie()
