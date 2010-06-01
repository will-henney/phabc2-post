from pyx import *
from math import pi

text.set(mode='latex')
unit.set(wscale=1.5, xscale=1.2)
figwidth = 10
figheight = figwidth/1.618
margin = 0.5


runtab = [
    ['95krum256', 1, 200, 20],
    ['95krumx256', 1, 600, 40],
    ]

skiphead = [0, 100]

datadir = '.'
istep = 1
tmax = 6			# in Myr

statsfile = [] 
vstatsfile = [] 
rstatsfile = [] 
dstatsfile = [] 
boxsizes = []
for runid, i1, i2, boxsize in runtab:
    statsfile += ['%s/%s-%4.4i-%4.4i-%4.4i.stats' % (datadir, runid, i1, i2, istep)]
    vstatsfile += ['%s/%s-%4.4i-%4.4i-%4.4i.vstats' % (datadir, runid, i1, i2, istep)]
    rstatsfile += ['%s/%s-%4.4i-%4.4i-%4.4i.rstats' % (datadir, runid, i1, i2, istep)]
    dstatsfile += ['%s/%s-%4.4i-%4.4i-%4.4i.dstats' % (datadir, runid, i1, i2, istep)]
    boxsizes += [boxsize]

mylines = graph.style.line()

## Graph 2 : ion frac, radius

c = canvas.canvas()

##

##
## Radius
##
d = []
for i in 0, 1:
    boxsize = boxsizes[i]
    d.append(graph.data.file(rstatsfile[i], x='Time/100', y='(boxsize/4.0)*rx2/3.085677582e18', context=locals(),
			     skiphead=skiphead[i], title=[r"$\left\langle R_\mathrm{ion}\right\rangle$", None][i]
			     ))
    d.append(graph.data.file(rstatsfile[i], x='Time/100', y='rif_min*boxsize/256', context=locals(), 
			     skiphead=skiphead[i], title=[r"$R_\mathrm{min}$", None][i]
			     ))
    d.append(graph.data.file(rstatsfile[i], x='Time/100', y='rif_max*boxsize/256', context=locals(), 
			     skiphead=skiphead[i], title=[r"$R_\mathrm{max}$", None][i]
			     ))
# d.append(graph.data.file(rstatsfile, x='Time/100', y='rmean_mass_i/3.086e18',
# 			 title=r"$\left\langle R\right\rangle_\mathrm{ion}$"
# 			 ))

g = graph.graphxy(width=figwidth, height=figheight, 
		  x=graph.axis.log(min=0.05, max=tmax, title='Time, Myr'), 
		  y=graph.axis.log(min=0.5, max=20,
				      title=r'Mean radius, parsec'),
		  key=graph.key.key(pos='tl', textattrs=[trafo.scale(0.7)],
				    hdist=0.3*unit.v_cm)
		  )

# homogeneous solution
# Where does this come from?
# R = R_0 (1 + 7 c_i t / 4 R_0)**4/7
# R_0 = (3 Q_H / 4 pi alpha n^2 )^1/3
#
# For the weak runs, we have Q_H = 5.e46 instead of 5.e48
# Also,  <T> is less: more like 8300 K (why?)
#
# Q_H = 5.e48 : R_0 = 0.539 pc with T=10^4 or 0.520 pc with T=9000
# Q_H = 5.e46 : R_0 = 0.116 pc with T=10^4 or 0.109 pc with T=8300
#
# rho c_i^2 = (1 + y_e) n k T => c_i = sqrt( (1 + y_e) k T / m )
#
# m = 1.3 mp, T = 8300 K, y_e = 1 => c_i = 1.027e6 cm/s (10.27 km/s)
#
# t_0 = 4 R_0 / 7 c_i = 54.4 (R_0 / pc) kyr
# 
# Actually, the above isn't quite right - best do it directly in python
# Run parameters
QH = 4.e46
#     Tmean = 8400.0
Tmean = 8900.0		# this gives the best fit to unif-weak-zerob256
n = 100.0
xi = 1.0
# physical constants
k = 1.3806503e-16
pc = 3.085677582e18
kyr = 1000.*3.15576e7
m = 1.3*1.67262158e-24
# This is taken from Garrelt's cgsconstants.f90
alpha = 2.59e-13*(Tmean/1.e4)**-0.7
# Find characteristic radius and time
R0 = (3*QH / (4*pi*alpha*n**2) )**(1./3.)
ci = ( (1.0+xi)*k*Tmean / m)**0.5
t0 = 4.*R0/(7.*ci) 
# Put in right units
R0 = R0 / pc
t0 = t0 / (1000*kyr)
print "R0 = %.2f pc, t0 = %.2f Myr" % (R0, t0)
f = graph.data.function(
    "y(x) = R0*(1.0 + x/t0)**(4./7.)", 
    context=locals(),
    title=r"$R_\mathrm{Str\ddot om}$"
)
g.plot(f, [graph.style.line([color.gray(0.8), style.linewidth.THIck])])
# plot the simulation lines last since they are thinner
threelines = attr.changelist([
	style.linestyle.solid, 
	style.linestyle.dashed, 
	style.linestyle.dotted, 
	])
g.plot(d, [graph.style.line([threelines])])

c.insert(g)

c.writePDFfile('radii_vs_t_' + 'krum')

##
## RMS and mean radial velocities
##
d = []

km = 1.e5

for i in 0, 1:
    boxsize = boxsizes[i]
    for Vel, title in [ 
			 ('Vrms_vol_i', 
			  r'$\left\langle v^2\right\rangle_\mathrm{ion}^{1/2}$'),
# 			 ('Vrms_vol_n', 
# 			  r'$\left\langle v^2\right\rangle_\mathrm{neut}^{1/2}$'), 
			 ]:
	d.append(graph.data.file(statsfile[i], x='Time/100', y=Vel+'/km', 
				 skiphead=skiphead[i], title=[title, None][i], context=locals()))

    for Vel, title in [ 
			 ('Vr_vol_i', r'$\left\langle v_r\right\rangle_\mathrm{ion}$'),
# 			 ('Vr_vol_n', r'$\left\langle v_r\right\rangle_\mathrm{neut}$'), 
			 ]:
	if i==1: title==None
	d.append(graph.data.file(vstatsfile[i], x='Time/100', y=Vel+'/km', 
				 skiphead=skiphead[i], title=[title, None][i], context=locals()))


g = graph.graphxy(width=figwidth, height=figheight,
		  x=graph.axis.linear(min=0, max=tmax, title='Time, Myr'), 
		  y=graph.axis.linear(min=0, max=3.0,
				      title=r'Mean gas velocities, km~s$^{-1}$'),
		  key=graph.key.key(pos='tr', textattrs=[trafo.scale(0.7)])
		  )

# homogeneous solution
# 
# V = (3/8) c_i (1 + 7 c_i t / 4 R_0)**-3/7
f = graph.data.function(
    "y(x) = (3./8.)*11.6*(1.0 + x/t0)**(-3./7.)",
    context=locals(),
    title=r'$\left\langle v_r\right\rangle_\mathrm{Str\ddot om}$')
g.plot(f, [graph.style.line([color.gray(0.8), style.linewidth.THIck])])

twolines = attr.changelist([
	style.linestyle.solid, 
	style.linestyle.dashed, 
	])
g.plot(d, [graph.style.line([twolines])])


g.writePDFfile('velocities_vs_t_' + 'krum')



