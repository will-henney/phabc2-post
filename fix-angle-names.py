import shutil
pattern = "t50m32-zerob-e01map-rot%+3.3i%+3.3i-%s0234.fits"
for i in range(360, 720, 5):
    j = i % 360
    for line in ['Halpha', 'N26584', 'O35007']:
	shutil.copy2(pattern % (i, i, line), pattern % (j, j, line))
