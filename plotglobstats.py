from pyx import *

text.set(mode='latex')
unit.set(wscale=1.5, xscale=1.2)


figwidth = 10
figheight = figwidth/1.618
margin = 0.5

datadir = '../../C2Ray/runs/globule/out/'
filetemplate = datadir + '/glob%s-%sx-globstats-table.dat'

runids = ['S90H', 'S90L', 'S00L', 'S45L', 'W90L', 'Z00L'] 
irun = {			
    # the old IDs
    'S90H': '12',
    'S90L': '12',
    'S00L': '11',
    'S45L': '05',
    'W90L': '14',
    'Z00L': '13',
    }

style = {			
    'S90H': graph.style.symbol(),
    'S90L': graph.style.line([style.linewidth.Thick]),
    'S00L': graph.style.line([style.linewidth.Thick, style.linestyle.dashed]),
    'S45L': graph.style.line([style.linewidth.Thick, style.linestyle.dotted]),
    'W90L': graph.style.line([style.linewidth.thin]),
    'Z00L': graph.style.line([style.linewidth.thin, style.linestyle.dashdotted]),
    }

vars = [
    # The individual variables we want to plot
    # [var, title, min, max, key position]
    ['Mn', r'Neutral mass, \(M_\odot\)', 0, 17, 'tr'],
    ['Mi', r'Ionized mass, \(M_\odot\)', 0, 17, 'tl'],
    ['Vnx', r'Neutral \(\langle V_x \rangle\), km s\(^{-1}\)', 0, 20, 'tl'],
    ['Vny', r'Neutral \(\langle V_y \rangle\), km s\(^{-1}\)', 0, 20, 'tl'],
    ['Bnx', r'Neutral \(\langle| B_x |\rangle\), \(\mu\)G', 0, 250, 'bl'],
    ['Bny', r'Neutral \(\langle| B_y |\rangle\), \(\mu\)G', 0, 250, 'bl'],
    ['Be_n', r'\(\langle\beta\rangle/\beta_0\), neutral gas', 0, 75, 'tl'],
    ['Be_i', r'\(\langle\beta\rangle/\beta_0\), ionized gas', 0, 75, 'bl'],
    ['Vxsh', r'Shocked shell \(\langle V_x \rangle\), km s\(^{-1}\)', 0, 20, 'tl'],
    ['Xif', r'Ionization front position \(x_\mathrm{if}\), pc', 0, 2, 'tl'],
    ['Xcom', r'Mean globule position \(\langle x\rangle\), pc', 0, 2, 'tl'],
    ]


for var, title, varmin, varmax, kpos in vars:
    g = graph.graphxy(width=figwidth, height=figheight, 
		      x=graph.axis.linear(title='Time, 1000~yr'), 
		      y=graph.axis.linear(min=varmin, max=varmax,
					  title=title),
		      key=graph.key.key(pos=kpos, textattrs=[trafo.scale(0.7)])
		      )
    for run in runids:
	if run.endswith('L'):
	    nz = '128'
	    every = 1
	    skiphead = 0
	    skiptail = 0
	else:
	    nz = '256'
	    every = 5
	    skiphead = 3
	    skiptail = 5
	datafile = filetemplate % (irun[run], nz)
	if var.startswith('Be'):
	    if run.startswith('W'):
		beta0 = 0.1
	    elif run.startswith('S'):
		beta0 = 0.01
	    else:
		continue 		# skip this one
	    d = graph.data.file(datafile, every=every, skiphead=skiphead, skiptail=skiptail, x='Time', y='%s/beta0' % var, title=run, context=locals())
	else:
	    d = graph.data.file(datafile, every=every, skiphead=skiphead, skiptail=skiptail, x='Time', y=var, title=run)

	g.plot(d, [style[run]])
    g.writePDFfile('glob-%s-stats.pdf' % (var))


