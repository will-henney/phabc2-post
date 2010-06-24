from movierot import Movie
import sys, os

def brightmax(i): 
    """
    Smoothly varying brightness with time so that the movie looks good

    This will have to be changed for the krum runs
    """
    return brightmax0/(1.0+(float(i)/100))


# Parse command line arguments
execname = os.path.split(sys.argv[0])[-1]
try: 
    runid = sys.argv[1]
    it1 = int(sys.argv[2])
    it2 = int(sys.argv[3])
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME1 ITIME2" % execname
    exit

if runid.find('krumx') >= 0:
    Movie.boxsize = 40.0
    brightmax0 = 5.e2
else:
    Movie.boxsize = 20.0
    brightmax0 = 4.e3

for itime in range(it1, it2+1):
    # Abuse the movie class to make one frame at a time
    movie = Movie(runid=runid, movieid="evo", 
		  time=itime, frame=itime, nframes=1, 
		  verbose=1, force=1)
    movie.brightmax = brightmax(itime)
    movie.camera.set_steps(0.0, 0.0) # don't move the camera!
    movie.bandscales = [0.4, 1.0, 0.05]
    movie.makeimages()

