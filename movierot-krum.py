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
    itime = int(sys.argv[2])
except IndexError, ValueError:
    print "Usage: %s RUNID ITIME" % execname
    exit

if runid.find('krumx') >= 0:
    Movie.boxsize = 40.0
    brightmax0 = 5.e2
else:
    Movie.boxsize = 20.0
    brightmax0 = 4.e3

th0, ph0 = 0.0, 0.0
movie = Movie(runid=runid, 
	      movieid="tumble%4.4i%+3.3i%+3.3i" % (itime, th0, ph0), 
	      time=itime, nframes=72, verbose=1, force=1)
movie.brightmax = brightmax(itime)
movie.camera.set_angles(th0, ph0)
movie.camera.set_steps(5.0, 5.0) # 72 frames for full revolution
movie.makemovie()
