from movierot import Movie
import sys, os

def brightmax(i): 
    """
    Smoothly varying brightness with time so that the movie looks good

    Changed 21 Jan 2008 - normalized to give 1.e7 @ 234,000 yrs
    for consistency with the tumble video

    TODO 21 Jan 2008 - Needs to be changed for the weak star 
    """
    return 1.7e7*(1.0+(float(i)/30)**2)/(1.0+(float(i)/50)**3)


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
movie = Movie(runid=runid, movieid="evo", 
	      time=it2, frame=1, nframes=it2-it1+1, 
	      verbose=1, force=-1) # force=-1 means don't make anything but movie
movie.makemovie()

