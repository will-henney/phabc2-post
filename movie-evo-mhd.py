from movierot import Movie
import sys, os

def brightmax(i): 
    """
    Smoothly varying brightness with time so that the movie looks good
    """
    return brightmax0*(1.0+(float(i)/30)**2)/(1.0+(float(i)/50)**3)


# Parse command line arguments
execname = os.path.split(sys.argv[0])[-1]
try: 
    runid = sys.argv[1]
    it1 = int(sys.argv[2])
    it2 = int(sys.argv[3])
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME1 ITIME2" % execname
    exit

Movie.boxsize = 4.0
if runid.find('weak') >= 0:
    brightmax0 = 5.e5
else:
    brightmax0 = 1.7e7

for itime in range(it1, it2+1):
    # Abuse the movie class to make one frame at a time
    movie = Movie(runid=runid, movieid="evo", 
		  time=itime, frame=itime, nframes=1, 
		  verbose=1, force=1)
    movie.brightmax = brightmax(itime)
    movie.camera.set_steps(0.0, 0.0) # don't move the camera!
    movie.makeimages()

