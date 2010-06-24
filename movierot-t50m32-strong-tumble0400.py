from movierot import Movie

def brightmax(i): 
    """
    Smoothly varying brightness with time so that the movie looks good
    """
    return brightmax0*(1.0+(float(i)/30)**2)/(1.0+(float(i)/50)**3)

brightmax0 = 2.e7
itime = 400

n = 9
dth = 360.0/n

movie = Movie(runid="t50m32dm0", movieid="tumble%4.4i" % itime, 
	      time=itime, nframes=n, verbose=1, force=1) # only recreate images, not maps
movie.brightmax = brightmax(itime)
movie.camera.set_steps(dth, dth)
movie.makemovie()


movie = Movie(runid="t50m32-zerob-e01", movieid="tumble%4.4i" % itime, 
	      time=itime, nframes=n, verbose=1, force=1) # only recreate images, not maps
movie.brightmax = brightmax(itime)
movie.camera.set_angles(360.0, 360.0)
movie.camera.set_steps(dth, dth)
movie.makemovie()
