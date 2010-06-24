from movierot import Movie

def brightmax(i): 
    """
    Smoothly varying brightness with time so that the movie looks good
    """
    return 2.6e7*(1.0+(float(i)/30)**2)/(1.0+(float(i)/50)**3)

it1, it2 = 100, 105
Movie.boxsize = 4.0
for itime in range(it1, it2+1):
    # Abuse the movie class to make one frame at a time
    movie = Movie(runid="t50m32e01", movieid="evo-test", 
		  time=itime, frame=itime, nframes=1, 
		  verbose=1, force=0)
    movie.makeimages()

movie = Movie(runid="t50m32e01", movieid="evo-test", 
	      time=it1, frame=it1, nframes=it2-it1+1, 
	      verbose=1, force=0) 
movie.makemovie()
