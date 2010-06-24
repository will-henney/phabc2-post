import sys, os
from movierot import Movie

brightmax0 = 2.e8

# Parse command line arguments
execname = os.path.split(sys.argv[0])[-1]
try: 
    runid = sys.argv[1]
    itime = int(sys.argv[2])
    brightmax0 = float(sys.argv[3])
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME" % execname
    exit

movie = Movie(runid=runid, movieid="tumble-%4.4i" % itime, 
	      time=itime, nframes=72, verbose=1, force=0)
movie.brightmax = brightmax0
movie.boxsize = 4.0
movie.camera.set_steps(5.0, 5.0)
movie.makemovie()
