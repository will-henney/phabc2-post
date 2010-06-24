from movierot import Movie

def brightmax(i): 
    """
    Smoothly varying brightness with time so that the movie looks good
    """
    return brightmax0*(1.0+(float(i)/30)**2)/(1.0+(float(i)/50)**3)

brightmax0 = 1.7e7
itime = 234

movie = Movie(runid="t50m32e01", movieid="tumble%4.4i" % itime, 
	      time=itime, nframes=72, verbose=1, force=2)
movie.brightmax = brightmax(itime)
movie.camera.set_steps(5.0, 5.0)
movie.makemovie()


movie = Movie(runid="t50m32-zerob-e01", movieid="tumble%4.4i" % itime, 
	      time=itime, nframes=72, verbose=1, force=2)
movie.brightmax = brightmax(itime)
movie.camera.set_angles(360.0, 360.0)
movie.camera.set_steps(5.0, 5.0)
movie.makemovie()
