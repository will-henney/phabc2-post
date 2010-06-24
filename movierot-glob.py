import sys, os
from movierot import Movie

def brightmax(i): 
    """
    Smoothly varying brightness with time so that the movie looks good
    """
    return brightmax0/(1.0+(float(i)/50)**2)

brightmax0 = 2.e8

# Parse command line arguments
execname = os.path.split(sys.argv[0])[-1]
try: 
    runid = sys.argv[1]
    itime = int(sys.argv[2])
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME" % execname
    exit

movie = Movie(runid=runid, movieid="tumble-%4.4i" % itime, 
	      time=itime, nframes=72, verbose=1, force=0)

movie.brightmax = brightmax(itime)
movie.camera.set_steps(5.0, 5.0)
movie.makemovie()
