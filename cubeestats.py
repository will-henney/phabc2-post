"""
Calculate energy statistics for simulation cubes

+ Thermal energy, magnetic energy, kinetic energy for each phase
+ The idea is to work out the energy transfer efficiencies from the radiation field

"""
# Imports from standard library
import os, sys, argparse
# Imports from third-party libraries
import numpy, pyfits
# Imports from my own libraries
import molfrac

execname = os.path.split(sys.argv[0])[-1]

# Parse command line arguments
parser = argparse.ArgumentParser(description="Generate table of energy statistics vs time")
parser.add_argument("--runid", "-r", type=str, default='Bstar-ep', 
                    help='ID for model run')
parser.add_argument("--i1", type=int, default=1, 
                    help='First save time')
parser.add_argument("--i2", type=int, default=400, 
                    help='Last save time')
parser.add_argument("--istep", type=int, default=1, 
                    help='Step between save times')
parser.add_argument("--boxsize", "-b", type=float, default=4.0, 
                    help="Size of simulation cube along x-axis in parsecs")
args = parser.parse_args()              # we can now use args.runid, args.itime


MP = 1.67262158e-24                       # Proton rest mass [cgs]
MU = 1.3                                  # Mean atomic mass
PC = 3.085677582e18                       # Parsec [cgs]
GAMMA = 5./3.

def load_fits(varid, itime):
    return pyfits.open('%s-%s%4.4i.fits' % (args.runid, varid, itime))['PRIMARY'].data

def energy_fmt(energies):
    "Format for printing the list of per-phase energies plus their sum"
    assert len(energies) == 3
    return "%9.2e %9.2e %9.2e" % tuple(energies) + " | %9.2e" % (sum(energies))

is_first_time = True 
for itime in range(args.i1, args.i2+1, args.istep):
    vx, vy, vz, xn, AV, dd = [load_fits(id, itime) for id in ["vx", "vy", "vz", "xn", "AV", "dd"]]
    bx, by, bz, pp = [load_fits(id, itime) for id in ["bx", "by", "bz", "pp"]]
    dn = dd / (MP*MU)      # number density

    # Do things that need to be done only once, but which depend on the cube size
    if is_first_time:     
        nz, ny, nx = dd.shape 
        z, y, x = numpy.ogrid[0:nz, 0:ny, 0:nx] # Cartesian grid coordinate arrays
        zc, yc, xc = [n/2 - 0.5 for n in [nz, ny, nx]] # Source assumed to be at grid center
        dVol = (args.boxsize*PC/nx)**3       # Assume cubic cells
        # rr is a grid of scalar radius from center (initially in units of cell size)
        rr = numpy.sqrt((x - xc)**2 + (y - yc)**2 + (z - zc)**2) 
        # [ux, uy, uz] is the unit radius vector from center of grid
        ux, uy, uz = (x - xc) / rr, (y - yc) / rr, (z - zc) / rr
        rr *= args.boxsize/nx        # Rescale radius to be in parsecs
        is_first_time = False        # Make sure not to do them again

    # Find weights for ionized, neutral, molecular phases
    xmol = molfrac.molfrac(AV)
    weights = [
        xn*xmol,                            # molecular
        xn*(1.0-xmol),                      # neutral
        1.0-xn                              # ionized
        ]

    # radial velocity is dot product of vector velocity with unit radius vector
    vr = vx*ux + vy*uy + vz*uz

    # Calculate energies
    ke, rke, the, me, mass, ee = [], [], [], [], [], []
    for w in weights: # each energy (ke, the, ...) has 3 values (1 for each phase)
        # kinetic energy
        ke.append(numpy.sum(w * 0.5*dd*(vx**2 + vy**2 + vz**2)) * dVol)
        # "outward radial" kinetic energy
        rke.append(numpy.sum(w * 0.5*dd*vr*numpy.abs(vr)) * dVol)
        # thermal energy
        the.append(numpy.sum(w * pp / (GAMMA - 1.0)) * dVol)
        # magnetic energy
        me.append(numpy.sum(w * 0.5*(bx**2 + by**2 + bz**2)) * dVol)
        # mass in each phase
        mass.append(numpy.sum(w * dd) * dVol)
        # total energy
        ee.append(ke[-1] + the[-1] + me[-1])


    print "="*60
    print "Time: ", itime
    print "Kinetic energies :", energy_fmt(ke)
    print "Radial KE        :", energy_fmt(rke)
    print "Thermal energies :", energy_fmt(the)
    print "Magnetic energies:", energy_fmt(me)
    print "."*57
    print "Total energies   :", energy_fmt(ee)
    print "Mass             :", energy_fmt(mass)

print "="*60
print "                 :", "%9s %9s %9s | %9s" % ("molec", "neut", "ion", "TOTAL")
print "="*60
