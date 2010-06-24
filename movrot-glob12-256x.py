from movierot import Movie

movie = Movie(runid="glob12-256x", movieid="tumblebright0033", 
	      time=33, nframes=72, verbose=1, force=0)
movie.brightmax = 1.5e8
movie.camera.set_steps(5.0, 5.0)
movie.makemovie()
