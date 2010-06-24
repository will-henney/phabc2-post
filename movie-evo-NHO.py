import sys, os
from movierot import Movie

def brightmax(i): 
    """
    Smoothly varying brightness with time so that the movie looks good
    """
    return brightmax0/(1.0+(float(i)/100)**2)


try: 
    runid = sys.argv[1]
    it1 = int(sys.argv[2])
    it2 = int(sys.argv[3])
except (IndexError, ValueError):
    print "Usage: %s RUNID ITIME1 ITIME2 [IDELTA THETA PHI BRIGHTMAX]" % execname
    exit

# Optional orientation arguments
try: 
    idelta = int(sys.argv[4])
except (IndexError, ValueError):
    idelta = 1

try: 
    th0 = float(sys.argv[5])
    ph0 = float(sys.argv[6])
except (IndexError, ValueError):
    th0, ph0 = 0.0, 0.0

# Optional brightness arguments
# set defaults
brightmax0 = 2.e8
bandscales = Movie.bandscales
try: 
    # 7th arg is scalar: use to set brightmax0
    brightmax0 = float(sys.argv[7])
except IndexError:
    pass
except ValueError:
    try:
        # 7th arg is triplet: use to set bandscales
        bandscales = [float(x) for x in sys.argv[7].split()]
    except ValueError:
        # otherwise, just use defaults
        pass

Movie.bandscales = bandscales
movie = Movie(runid=runid, movieid="evo%+3.3i%+3.3i" % (th0, ph0), 
	      time=it1, dtime=idelta, frame=it1/idelta, 
              nframes=(it2-it1)/idelta + 1, verbose=1, force=0)
movie.brightmaxfunc = brightmax
movie.camera.set_steps(0.0, 0.0)
movie.camera.set_angles(th0, ph0)
movie.boxsize = 4.0
movie.makemovie()
