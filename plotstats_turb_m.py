from pyx import *

text.set(mode='latex')
text.preamble(r"\usepackage[varg]{txfonts}")
unit.set(wscale=1.5, xscale=1.2)
figwidth = 10
figheight = figwidth/1.618
margin = 0.5


runtab = [
    # ['Ostar', 1, 92 ],
    ['Bstar-ep', 10, 1020 ],
    ['Bstar-HDep', 10, 1000 ],
    ]


datadir = '.'
istep = 10
tmax = 1200

for runid, i1, i2 in runtab:
    print "Creating stats graph for ", runid
    statsfile = '%s/%s-%4.4i-%4.4i-%4.4i.stats' % (datadir, runid, i1, i2, istep)
    vstatsfile = '%s/%s-%4.4i-%4.4i-%4.4i.vstats' % (datadir, runid, i1, i2, istep)
    rstatsfile = '%s/%s-%4.4i-%4.4i-%4.4i.rstats' % (datadir, runid, i1, i2, istep)
    dstatsfile = '%s/%s-%4.4i-%4.4i-%4.4i.dstats' % (datadir, runid, i1, i2, istep)

    mylines = graph.style.line()

    ## Graph 1 : density, clumping

    c = canvas.canvas()

    ##
    ## Mean density
    ##
    d = []

    for Dmean, title in [ 
	('Dmean_i', r'$\langle n\rangle_\mathrm{ion}$'),
	('Dmean_n', r'$\langle n\rangle_\mathrm{neut}$'), 
	('Dmean_tot', r'$\langle n\rangle_\mathrm{tot}$'), 
	]:
	d.append(graph.data.file(dstatsfile, x='Time', y=Dmean, title=title))

    g = graph.graphxy(width=figwidth, height=figheight, 
    # 		  x=graph.axis.linear(title='Time, 1000~yr'), 
		      x=graph.axis.linear(min=0, max=tmax, painter=graph.axis.painter.linked()),
		      y=graph.axis.logarithmic(min=5, max=5.e5,
					  title=r'Mean densities, cm$^{-3}$'),
		      key=graph.key.key(pos='tr', textattrs=[trafo.scale(0.7)])
		      )
    g.plot(d, [graph.style.line()])

    c.insert(g, [trafo.translate(0, figheight+margin)])

    ##
    ## Clumping 
    ##
    d = []
    for Clump, title in [ 
			  ('D2mean_i/Dmean_i**2', r'$C_\mathrm{ion}$'),
			  ('D2mean_n/Dmean_n**2', r'$C_\mathrm{neut}$'), 
			  ('D3mean_i**2/D2mean_i**3', r'$\varepsilon_\mathrm{ion}^{-1}$')]:
	d.append(graph.data.file(dstatsfile, x='Time', y=Clump, title=title))

	g = graph.graphxy(width=figwidth, height=figheight, 
			  x=graph.axis.linear(min=0, max=tmax, title='Time, 1000~yr'), 
			  y=graph.axis.logarithmic(title=r'Degree of clumping',
    # 	    r'Clumping factor, $C = \langle n^2 \rangle / \langle n \rangle^2$',
						   min=0.9, max=450
		),
		      key=graph.key.key(pos='tr', textattrs=[trafo.scale(0.7)])
		      )

    g.plot(d, [mylines])

    g.writePDFfile('clumping_vs_t_' + runid)

    c.insert(g)

    c.writePDFfile('stats1_vs_t_' + runid)


    ## Graph 2 : ion frac, radius

    c = canvas.canvas()

    ##
    ## Ion frac
    ##
    d = []

    for Frac, title in [ ('Ifrac_v2', r'$X_\mathrm{vol}$'), 
			 ('Ifrac_m', r'$X_\mathrm{mass}$') ]:
	d.append(graph.data.file(statsfile, x='Time', y=Frac, title=title))

    g = graph.graphxy(width=figwidth, height=figheight, 
    # 		  x=graph.axis.linear(title='Time, 1000~yr'), 
		      x=graph.axis.linear(min=0, max=tmax, painter=graph.axis.painter.linked()),
		      y=graph.axis.logarithmic(min=1.e-6, max=1, title=r'Total ionized fraction'),
		      key=graph.key.key(pos='tl', textattrs=[trafo.scale(0.7)])
		      )

    g.plot(d, [graph.style.line()])

    c.insert(g, [trafo.translate(0, figheight+margin)])


    ##

    ##
    ## Radius
    ##
    d = []
    d.append(graph.data.file(rstatsfile, x='Time', y='rx2/3.086e18', 
			     title=r"$\left\langle R_\mathrm{ion}\right\rangle$"
			     ))
    d.append(graph.data.file(rstatsfile, x='Time', y='rif_min*4.0/256', 
			     title=r"$R_\mathrm{min}$"
			     ))
    d.append(graph.data.file(rstatsfile, x='Time', y='rif_max*4.0/256', 
			     title=r"$R_\mathrm{max}$"
			     ))
    # d.append(graph.data.file(rstatsfile, x='Time', y='rmean_mass_i/3.086e18',
    # 			 title=r"$\left\langle R\right\rangle_\mathrm{ion}$"
    # 			 ))

    g = graph.graphxy(width=figwidth, height=figheight, 
		      x=graph.axis.linear(min=0, max=tmax, title='Time, 1000~yr'), 
		      y=graph.axis.linear(min=0, max=4,
					  title=r'Mean radius, parsec'),
		      key=graph.key.key(pos='tl', textattrs=[trafo.scale(0.7)],
					hdist=0.3*unit.v_cm)
		      )

    g.plot(d, [graph.style.line()])

    # homogeneous solution
    f = graph.data.function(
	"y(x) = 0.532*(1.0 + 0.0390*x)**(4./7.)", 
	title=r"$R_\mathrm{Str\ddot om}$"
    )
    g.plot(f, [graph.style.line([color.gray(0.8), style.linewidth.Thick])])

    c.insert(g)

    c.writePDFfile('stats2_vs_t_' + runid)

    ##
    ## RMS and mean radial velocities
    ##
    d = []

    km = 1.e5

    for Vel, title in [ 
			 ('Vrms_vol_i', 
			  r'$\left\langle v^2\right\rangle_\mathrm{ion}^{1/2}$'),
			 ('Vrms_vol_n', 
			  r'$\left\langle v^2\right\rangle_\mathrm{neut}^{1/2}$'), 
			 ]:
	d.append(graph.data.file(statsfile, x='Time', y=Vel+'/km', title=title, context=locals()))

    for Vel, title in [ 
			 ('Vr_vol_i', r'$\left\langle v_r\right\rangle_\mathrm{ion}$'),
			 ('Vr_vol_n', r'$\left\langle v_r\right\rangle_\mathrm{neut}$'), 
			 ]:
	d.append(graph.data.file(vstatsfile, x='Time', y=Vel+'/km', title=title, context=locals()))


    g = graph.graphxy(width=figwidth, height=figwidth,
		      x=graph.axis.linear(min=0, max=tmax, title='Time, 1000~yr'), 
		      y=graph.axis.linear(min=0, max=13,
					  title=r'Mean gas velocities, km~s$^{-1}$'),
		      key=graph.key.key(pos='tr', textattrs=[trafo.scale(0.7)])
		      )

    g.plot(d, [graph.style.line()])

    # homogeneous solution
    f = graph.data.function(
	"y(x) = (3./8.)*11.6*(1.0 + 0.0390*x)**(-3./7.)",
	title=r'$\left\langle v_r\right\rangle_\mathrm{Str\ddot om}$')
    g.plot(f, [graph.style.line([color.gray(0.8), style.linewidth.Thick])])

    g.writePDFfile('stats3_vs_t_' + runid)



