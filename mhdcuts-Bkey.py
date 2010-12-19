from PIL import Image
import colorsys
import numpy as np
nx, ny = 128, 128

def createHSVimage(hue, sat, val, alpha):
    rgba = []
    # convert hsv -> rgb
    for h, s, v, a in zip(hue.tolist(), sat.tolist(), val.tolist(), alpha.tolist()):
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        rgba.append( (int(255*r), int(255*g), int(255*b), int(255*a)) )
    # new image for color composite
    imrgb = Image.new('RGBA', (nx, ny))
    imrgb.putdata(rgba)
    imrgb = imrgb.transpose(Image.FLIP_TOP_BOTTOM)
    return imrgb

x = np.linspace(-1.2, 1.2, nx)[np.newaxis,:]
y = np.linspace(-1.2, 1.2, ny)[:,np.newaxis]

rr = np.sqrt(x**2 + y**2)
th = np.arctan2(y, x)

phase = 0.0
boost = 2.0			# exaggerate the field angle
phi = np.remainder(phase + boost*np.arctan2(x, y), 
                  2.0*np.pi)/(2.0*np.pi) # should be in range [0, 1]
# dip = np.where(rr <= 1.0, 1.0 - np.sqrt(1.0 - rr**2), 1.0)
dip = np.where(rr <= 1.0, rr**2, 1.0)
hue = phi

bgamma = 2.0
# hue is in-plane angle
hue = phi.flatten()
# saturation is out-of-plane angle
sat = dip.flatten()

# alpha channel to smoothly cut out the circle
sharpness = 20.0
alpha = 1.0/(1.0 + np.exp(sharpness*(rr-1.05))).flatten()

for field_strength in [0.5, 1.0]:
    # value is field strength
    val = np.array([field_strength]*nx*ny)
    imrgb = createHSVimage(hue, sat, val, alpha)
    imrgb.save("mhdcuts-Bkey-%2.2i.png" % (int(10*field_strength)))

# and another one, with saturation going down in the middle but value too
val = np.where(rr <= 1.0, 0.8 + 0.2*rr**2, 1.0).flatten()
imrgb = createHSVimage(hue, sat, val, alpha)
imrgb.save("mhdcuts-Bkey-sat.png") 

# now do another one, with full saturation but value proportional to radius
sat = np.ones(nx*ny, np.float)
val = np.where(rr <= 1.0, rr, 1.0).flatten()
imrgb = createHSVimage(hue, sat, val, alpha)
imrgb.save("mhdcuts-Bkey-val.png") 

# and finally one that combines everything....
peak = 0.7                              # fully saturated, value=1 appears here
val = np.where(rr >= peak,  (peak/rr)**2 , 1.0).flatten()
sat = np.where(rr <= peak, (rr/peak)**3, 1.0).flatten()
imrgb = createHSVimage(hue, sat, val, alpha)
imrgb.save("mhdcuts-Bkey-both.png") 
