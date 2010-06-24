import sys, os
from movierot import Movie

def brightmax(i): 
    """
    Smoothly varying brightness with time so that the movie looks good
    """
    return brightmax0/(1.0+(float(i)/50)**2)


# Parse command line arguments
execname = os.path.split(sys.argv[0])[-1]
try: 
    runid = sys.argv[1]
    it1 = int(sys.argv[2])
    it2 = int(sys.argv[3])
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME1 ITIME2 [THETA PHI BRIGHTMAX]" % execname
    exit

# Optional orientation arguments
try: 
    th0 = float(sys.argv[4])
    ph0 = float(sys.argv[5])
except IndexError, ValueError:
    th0, ph0 = 0.0, 0.0
# Optional brightness arguments
try: 
    brightmax0 = float(sys.argv[6])
except IndexError, ValueError:
    brightmax0 = 2.e8


movie = Movie(runid=runid, movieid="evo%+3.3i%+3.3i" % (th0, ph0), 
	      time=it1, frame=it1, nframes=it2-it1+1, verbose=1, force=0)
movie.brightmaxfunc = brightmax
movie.camera.set_steps(0.0, 0.0)
movie.camera.set_angles(th0, ph0)
movie.dtime = 1

movie.makemovie()
