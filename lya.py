"""
Create all the datacubes required for Lyman-alpha transport

+ Temperature, T
+ Velocity field, *v*
+ Density of neutral hydrogen, n_HI
+ Density of ionized hydrogen, n_HII
+ Metallicity, Z
+ Lya emissivity, L_Lya
+ Continuum emissivity around the Lya line, L_cont

As requested by Peter Laursen <pela@dark-cosmology.dk>
Via Garrelt Mellema <garrelt@astro.su.se>

"""
import argparse
import os
import numpy
from numpy import pi
import pyfits
import molfrac

# WJH : this doesn't work because fortran-namelist/namelist.py is broken
# import namelist
# class Simulation(object):
#     attrs = ["qh0", "tstar"]
#     def __init__(self, runid):
#         nl = namelist.Namelist("../in/%s.pars" % (runid))

class Run(object):
    pass
        
# Parameters from Bstar-ep.pars
# https://github.com/deprecated/phabc2-runpars/blob/master/pdrturb/in/Bstar-ep.pars
Run.xmax =  1.2344e19   
Run.ymax =  1.2344e19   
Run.zmax =  1.2344e19   
Run.time_save_interval =  3.15576e10 
Run.nxg = 256
Run.nyg = 256
Run.nzg = 256  
Run.i0 = 128
Run.j0 = 128
Run.k0 = 128

bstar = Run()
bstar.qh0 = 4.e46
bstar.tstar = 3.3e+04
bstar.fuv_euv_ratio = 62.5
bstar.nsave_total_run = 2000    

ostar = Run()
ostar.qh0 = 5.e48
ostar.tstar = 3.75e+04,
ostar.fuv_euv_ratio = 6.25
ostar.nsave_total_run = 400    

MP = 1.67262158e-24
MU = 1.3
PC = 3.085677582e18

def lya_emissivity(ni, T):
    """
    True emissivity of Lyman-alpha in photons/cm^3/s/sr

    Only includes recombination contribution. Does not include
    collisional excitation or scattering or fluorescence.

    Calculated in Case A approximation
    """
    # Data from Osterbrock & Ferland 
    alpha_eff = 2.04e-14 * 32.7 * (1216.0/4863.0) * (1.e4/T)
    # this gives 1.67e-13 @ 1.e4 K
    return alpha_eff * ni*ni / (4*pi)

def fuv_field(av, r):
    """
    Mean intensity of FUV field around 1216 Ang

    This is actually the total photon intensity in the range 912 to 2000 Ang

    Units are photons/cm^2/s/sr

    Input arguments: 

    av - V-band extinction from star to local point (mag)
    r - distance from star to local point (cm)
    """
    return numpy.exp(-1.9*av)/(4.0*pi*r*r)

def load_fits(id):
    return pyfits.open(
        '%s-%s%4.4i.fits' % (args.runid, id, args.itime)
        )['PRIMARY'].data

def write_fits(fitsfile, data):
    hdu = pyfits.PrimaryHDU(data)
    hdu.writeto(savedir + "/" + fitsfile + ".fits", clobber=True)

parser = argparse.ArgumentParser(description="Create cubes for Lyman alpha")
parser.add_argument("--runid", "-r", type=str, default='Bstar-ep', 
                    help='ID for model run')
parser.add_argument("--itime", "-i", type=int, default=1, 
                    help='Integer save time counter')
args = parser.parse_args() 

if args.runid.startswith("Ostar"):
    star = ostar
elif args.runid.startswith("Bstar"):
    star = bstar
else:
    raise NotImplementedError("Only Ostar and Bstar are implemented, sorry!")

savedir = "lya-%s-%4.4i" % (args.runid, args.itime)

if not os.path.isdir(savedir):
    os.mkdir(savedir)

vx, vy, vz, xn, av, dd, te = [load_fits(id) for id in 
                              ["vx", "vy", "vz", "xn", "AV", "dd", "te"]]

x, y, z, = numpy.ogrid[0:star.nxg, 0:star.nyg, 0:star.nzg]       # grid coordinates
# grid of scalar radius from center (in units of cell size)
rr = numpy.sqrt((x - star.i0)**2 + (y - star.j0)**2 + (z - star.k0)**2)
rr[rr==0] = 1.0
rr *= star.xmax/star.nxg                # convert to cm

mm = molfrac.molfrac(av)
ni = dd*(1.-xn) / MP / MU
nn = dd*xn*(1.-mm) / MP / MU
nm = dd*xn*mm / MP / MU
ly = lya_emissivity(ni, te)
jc = fuv_field(av, rr)
print jc.max(), jc.min()
jc *= star.qh0 * star.fuv_euv_ratio / (4*pi)
print jc.max(), jc.min()

for id, data in [ 
    ["temperature", te],
    ["vel_x", vx],
    ["vel_y", vy],
    ["vel_z", vz],
    ["n_hi", nn],
    ["n_hii", ni],
    ["n_h2_x_2", nm],
    ["lya_emissivity", ly],
    ["continuum_intensity", jc],
    ]:
    write_fits(id, data)



