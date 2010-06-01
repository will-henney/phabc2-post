"""
First stab at redoing the appendix plots using PyX

09 Apr 2009
"""

import heatcool
import numpy as N

import pyx
pyx.text.set(mode="latex")

import os

####
# 
# General setup
# 
####

import cloudy
cloudy.datadir = os.path.join(os.getenv("HOME"), "Work", "JaneCool", "out")
cloudy_file_template = "pdr-new-h2x-D%.1f-n%.1f-%s.ovr"

figwidth = 10
margin = 0.2

gridstyles = [pyx.style.linewidth.thin, pyx.color.transparency(0.8)]
mypainter = pyx.graph.axis.painter.regular(gridattrs=gridstyles)
linkpainter = pyx.graph.axis.painter.linked(gridattrs=gridstyles)

basename = __file__.split('.')[0]

####
# 
# Cooling graphs
# 
####

Tmin, Tmax = 10.0, 1.5e4
Coolmin, Coolmax = 1.e-34, 1.e-23

gleft = pyx.graph.graphxy(width=figwidth, 
			  x=pyx.graph.axis.log(title=r"Gas temperature, K", 
					       min=Tmin, max=Tmax, painter=mypainter), 
			  y=pyx.graph.axis.log(title=r"Cooling coefficients, \(L/n^2\), erg cm\(^3\) s\(^{-1}\)", 
					       min=Coolmin, max=Coolmax, painter=mypainter),
			  key=pyx.graph.key.key(pos="br", 
						keyattrs=[pyx.deco.filled, pyx.color.rgb.white], 
						textattrs=[pyx.trafo.scale(0.7)]))

gright = pyx.graph.graphxy(width=figwidth, 
			  x=pyx.graph.axis.log(title=r"Gas temperature, K", 
					       min=Tmin, max=Tmax, painter=mypainter), 
			  y=pyx.graph.axis.linkedaxis(gleft.axes["y"], painter=linkpainter),
			  key=pyx.graph.key.key(pos="br", 
						keyattrs=[pyx.deco.filled, pyx.color.rgb.white], 
						textattrs=[pyx.trafo.scale(0.7)]))

cloudy_linestyles = [pyx.style.linewidth.thick, 
		      pyx.style.linestyle.solid, 
		      pyx.style.linejoin.round, 
		      pyx.color.transparency(0.3)]

fit_linestyles = [pyx.style.linewidth.THick, 
		  pyx.style.linestyle.solid, 
		  pyx.style.linejoin.round, 
		  pyx.color.transparency(0.3)]


density_colorstyles = {
    2.0: [pyx.color.cmyk.Dandelion],
    4.0: [pyx.color.cmyk.SeaGreen],
    6.0: [pyx.color.cmyk.RedViolet],
    }

for logden in [2.0, 4.0, 6.0]:
    first_in_batch = True
    cloudylines = []
    for D in [0.1, 0.3, 0.5, 1.0]:
	for radstring in ["XO", "O"]:
	    cloudy_file = cloudy_file_template % (D, logden, radstring)
	    dataset = cloudy.datasetfile(cloudy_file)
	    Temperature = 10**dataset.grabcolumn('Te')
	    Cool = 10**(dataset.grabcolumn('Htot') - 2.0*dataset.grabcolumn('hden'))
	    if first_in_batch:
		title = r"Cloudy models: \(n = 10^{%d}~\mathrm{cm}^{-3}\)" % (logden)
		first_in_batch = False
	    else:
		title = None
	    cloudylines.append(
		pyx.graph.data.values(x=Temperature, y=Cool, title=title)
		)
    # plot all the same density in the same style
    gleft.plot( cloudylines, 
		[pyx.graph.style.line(cloudy_linestyles + density_colorstyles[logden])] 
		)

    Temperature = N.logspace(1.2, 4.0, num=50)
    Cool = heatcool.cool(10**logden, Temperature)
    fitline = pyx.graph.data.values(x=Temperature, y=Cool, title=r"Fit: \(n = 10^{%d}~\mathrm{cm}^{-3}\)" % (logden))
    gright.plot(fitline,  [ 
	    pyx.graph.style.line(fit_linestyles + density_colorstyles[logden]),
	    ] )


# Finally, plot the KI2002 curve
Temperature = N.logspace(1.2, 4.0, num=50)
Cool = heatcool.coolKI(Temperature)
fitline = pyx.graph.data.values(x=Temperature, y=Cool, title=r"Koyama \& Inutsuka (2002)")
gright.plot(fitline,  [ 
	pyx.graph.style.line(fit_linestyles + [pyx.color.cmyk.Lavender, pyx.style.linestyle.dashed]),
	] )

c = pyx.canvas.canvas()
c.insert(gleft)
c.insert(gright, [pyx.trafo.translate(figwidth + margin, 0)])

c.writePDFfile("%s-cool" % (basename))

####
# 
# Heating graphs
# 
####

Avmin, Avmax = 0.0+1.e-5, 20.0
Heatmin, Heatmax = 1.e-28, 1.e-20

heatcool.cosmicrayfactor = 0.5

D_styles = {
    0.1: [pyx.style.linewidth.THick],
    0.5: [pyx.style.linewidth.normal],
    }

X_styles = {
    "XO": [pyx.style.linestyle.solid],
    "O": [pyx.style.linestyle.dashed],
    }

cloudy_linestyles = [pyx.style.linewidth.normal, 
		      pyx.style.linestyle.solid, 
		      pyx.style.linejoin.round, 
		      pyx.color.transparency(0.3)]

fit_linestyles = [pyx.style.linewidth.normal, 
		  pyx.style.linestyle.solid, 
		  pyx.style.linejoin.round, 
		  pyx.color.transparency(0.3)]

gleft = pyx.graph.graphxy(width=figwidth, 
			  x=pyx.graph.axis.lin(title=r"Visual extinction, \(A_V\), magnitudes", 
					       min=Avmin, max=Avmax, painter=mypainter), 
			  y=pyx.graph.axis.log(title=r"Heating per particle, \(H/n\), erg s\(^{-1}\)", 
					       min=Heatmin, max=Heatmax, painter=mypainter),
			  key=pyx.graph.key.key(pos="tr", 
						keyattrs=[pyx.deco.filled, pyx.color.rgb.white], 
						textattrs=[pyx.trafo.scale(0.7)]))

gright = pyx.graph.graphxy(width=figwidth, 
			  x=pyx.graph.axis.lin(title=r"Visual extinction, \(A_V\), magnitudes", 
					       min=Avmin, max=Avmax, painter=mypainter), 
			  y=pyx.graph.axis.linkedaxis(gleft.axes["y"], painter=linkpainter),
			  key=pyx.graph.key.key(pos="tr", 
						keyattrs=[pyx.deco.filled, pyx.color.rgb.white], 
						textattrs=[pyx.trafo.scale(0.7)]))

for D in [0.1, 0.5]:
    for radstring in ["XO", "O"]:
	for logden in [2.0, 4.0, 6.0]:
	    if "X" in radstring:
		title = r"\(r = %.1f~\mathrm{pc}\), \(n = 10^{%d}~\mathrm{cm}^{-3}\)" % (D, logden)
	    else:
		title = None
	    cloudy_file = cloudy_file_template % (D, logden, radstring)
	    dataset = cloudy.datasetfile(cloudy_file)
	    Av = dataset.grabcolumn('AV(point)')
	    Heat = 10**(dataset.grabcolumn('Htot') - dataset.grabcolumn('hden'))
	    cloudyline = pyx.graph.data.values(x=Av, y=Heat, title=title)
	    gleft.plot( cloudyline, 
			[pyx.graph.style.line(cloudy_linestyles 
					      + D_styles[D] 
					      + X_styles[radstring]
					      + density_colorstyles[logden]
					      )] 
			)
	    Av = N.linspace(0.0, 20.0, num=50)
	    Heat = heatcool.heat(10**logden, Av, D, radiation=radstring)
	    fitline = pyx.graph.data.values(x=Av, y=Heat, title=title)
	    gright.plot(fitline,  
		    [pyx.graph.style.line(fit_linestyles 
					  + D_styles[D] 
					  + X_styles[radstring]
					  + density_colorstyles[logden]
					  )] 
		    )

c = pyx.canvas.canvas()
c.insert(gleft)
c.insert(gright, [pyx.trafo.translate(figwidth + margin, 0)])

c.writePDFfile("%s-heat" % (basename))
