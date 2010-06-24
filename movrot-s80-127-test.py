from movierot import Movie

movie = Movie(runid="s80-127", movieid="tumble0060", 
	      time=60, nframes=72, verbose=1, force=0)
movie.brightmax = 1.5e8
movie.camera.set_steps(5.0, 5.0)
movie.makemovie()
