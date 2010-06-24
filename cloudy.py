#! /usr/bin/env python
##
## cloudy.py: (c) Will Henney  2004
##
## Classes and functions for parsing cloudy data file and dealing with the
## resulting data structures

import numpy as numarray # used by heating/cooling parser
import re, os

# module variables (may be overwridden by user)
datadir = './'

##
## Some generic utility functions
##
def maybefloat( x ):
    "Do float(x) if possible. If not, just x"
    try:
        return float(x)
    except:
        return x

##
## New paste algorithm 05 Mar 2004 - cleaner, no dependency on
## Numeric, plus works with all data types - still missing error
## checks though
##
def paste ( tab_a, tab_b ) :
    """Paste together the columns of a 2D table
    - this is modeled on the Unix paste command
    - each table must have same number of rows"""
    tab_c = []
    for irow in range(len(tab_a)) :
        tab_c.append( tab_a[irow] + tab_b[irow] )
    return tab_c


##
## Basic class idea is copied from similar structures in pyx/data.py
## (but vastly simplified)
##

class dataset:
    "Minimal Cloudy dataset class"

    titles = []
    """column titles (list of strings)"""

    points = []
    """column/row/iteration data
    - a list of iterations
    - each iteration is a list of rows (each row represents a spatial point)
    - each row is a list of floats, corresponding to each member of titles

    Eventually I may define methods to extract iterations and the like, but
    for the time being you have to use the normal list notation. E.g.,
    last_iteration = dataset.data[-1]"""

    iter_separator = '\n\n' # 2 blank lines means 3 newlines
    row_separator = '\n'
    field_separator = '\t'
    comment_character = '#' 
    def __init__(self, points=[], titles=[]):
        """initializes a dataset instance
        - just here for debugging really
        - in real use one would use the datasetfile child class"""
        self.points = points
        self.titles = titles

    def writefile(self, file, itnumber=-1):
        try:
            f = open(os.path.join(datadir, file), "w")
        except:
            print "Error opening cloudy data file"
	    raise IOError
        headerline = self.comment_character \
	    + self.field_separator.join([title for title in self.titles]) \
	    + '\n'
        f.write(headerline)
        datalines = self.row_separator.join( 
	    [self.field_separator.join(["%7e" % x for x in line]) 
	     for line in self.points[itnumber] ] )
        f.write(datalines)
        f.close()

class dumbdatafile(dataset) :
    """This is a class that just has a list of iterations, each of which
    can be written to a separate file"""

    iterations = []
    headerstring = ""

    def __init__(self, file):
        "Read data from a file"
        try:
            f = open(os.path.join(datadir, file), "r")
        except:
            print "Error opening cloudy data file"
	    raise IOError
        # get all the lines from the file as a single string
        lines = f.read()
        self.iterations = lines.split(self.iter_separator)
        self.headerstring = self.iterations[0].split(self.row_separator)[0]

    def writefile(self, file, itnumber=-1):
        try:
            f = open(datadir + file, "w")
        except:
            print "Error opening cloudy data file"
	    raise IOError
        f.write(self.headerstring)
        f.write(self.iterations[itnumber])
        f.close()

class datasetfile(dumbdatafile):
    "A Cloudy dataset class that is initialized from a file"

    def __init__(self, file):
        dumbdatafile.__init__(self,file)
        self.points = []
        for iteration in self.iterations:
            rows = iteration.split(self.row_separator)
            iteration = [] # nuke this iteration so we can fill it again
            for row in rows:
                try: 
                    if row.startswith('#'):
#                     if row.startswith('#') and row.find('depth') >= 0:
#                         # Take the column headers from the relevant
#                         # comment line after first removing the
#                         # `#'. Note that this assumes that the first
#                         # column is always going to be depth.
                        self.titles = row[1:].split(self.field_separator)
                    elif row and not row.startswith('#'):
                        # All other non-empty, non-comment rows are
                        # appended to the current iteration. Each field
                        # is converted from string to float if
                        # possible
                        iteration.append(  
			    [maybefloat(field) for field in 
			     row.split(self.field_separator)])
                except:
                    # Failed for some reason - probably an empty row
                    # Just ignore it
                    continue
            # now that this iteration is sorted, we can stuff it into points
            self.points.append( iteration )

    def grabcolumn( self, columnname, iter=-1, firstrow=0, nrows=0) :
	""" Grabs a whole column from a cloudy.datasetfile instance as a
	numarray vector.  By default, the last iteration is used.  """
	if nrows:
	    lastrow = firstrow + nrows -1
	else :
	    lastrow = -1
	try :
	    column = []
	    for row in self.points[iter][firstrow:lastrow]:
		value = row[self.titles.index(columnname)]
		if isinstance(value, float): 
		    column.append(value)
	    return numarray.array(column)
	except IndexError:
	    print "Warning: Failed to grab column!"
	    return numarray.zeros(len(t),numarray.Float)

    def grabtextcolumn( self, columnname, iter=-1, firstrow=0, nrows=0) :
	""" Grabs a whole column from a cloudy.datasetfile instance as a
	numarray vector.  By default, the last iteration is used.  """
	if nrows:
	    lastrow = firstrow + nrows -1
	else :
	    lastrow = -1
	try :
	    column = []
	    for row in self.points[iter][firstrow:lastrow]:
		column.append(row[self.titles.index(columnname)])
	    return column
	except IndexError:
	    print "Warning: Failed to grab column!"
	    return ''*len(self.points[iter][firstrow:lastrow])
	

def cleanup_dr_file ( filename ) :
    """Cleans up a .dr file written by cloudy in order to:

1. Insert double new lines between each iteration
2. Remove the garbage lines such as ' >>>> A temperature failure occured.'
3. Remove the first few lines, with cell 0
"""
    try:
        f = open(datadir + filename, "r")
    except:
        print "Error opening cloudy data file"
    drlines = f.readlines()
    f.close()
    f = open(datadir + filename+".new", "w")
    f.write(drlines[0]) # write out comment line
    for line in drlines[1:]:
        if line.startswith("1\t"):
            f.write("\n\n")
        if re.match("^[1-9]",line):
            f.write(line)
    f.close()    



def findgradient( x, y, n=1 ):
    dydx = numarray.zeros(len(x),numarray.Float)
    # loop over zones that are far enough from endpoints
    nbad = 0
    nxybad = 0
    nx2bad = 0
    for i in range(n,len(x)-n):
        # x is w.r.t. mean x
        xbar = numarray.sum(x[i-n:i+n+1])/(2.*n+1.)
        xwin = (x[i-n:i+n+1]-xbar)/xbar
        # y is w.r.t. mean y
        ybar = numarray.sum(y[i-n:i+n+1])/(2.*n+1.)
        ywin = (y[i-n:i+n+1]-ybar)/ybar
        sumxy = numarray.sum(xwin*ywin)
        sumx2 = numarray.sum(xwin**2)
        try:
            dydx[i] = sumxy/sumx2
        except ZeroDivisionError, FloatingPointError:
            nbad+=1
            if sumx2==0.0:
                nx2bad+=1
            if sumxy==0.0:
                nxybad+=1
            dydx[i] = 0.0
        dydx[i] *= ybar/xbar # put back in normalization constants
    if nbad:
        print 'Warning: failed to calculate dy/dx in %i zones out ot %i' \
	    % (nbad, len(x))
        print '         sumxy = 0.0 in %i zones' % nxybad
        print '         sumx2 = 0.0 in %i zones' % nx2bad

    # fill in end zones
    dydx[:n] = dydx[n]
    dydx[-n:] = dydx[-n-1]
    return dydx
    
class enthalpyset(dataset):
    "Dataset that contains the specific enthalpy and its derivative"

    # assume standard abundances
    amass = 1.3*1.67e-24 # mean atomic mass in grams

    def __init__(self, prefix, n=1) :
        """Reads both the pressure file and the drad file, which are necessary
        for calculating all the enthalpy quantities"""

        self.points = []
        
        # grab the pressure and dr data
        print "Reading pressure file"
        self.pdatafile = datasetfile(prefix+".pre")
        print "Reading dr file"
        self.drdatafile = datasetfile(prefix+".dr")
        print "Reading overview file"
        self.odatafile = datasetfile(prefix+".ovr")

        niter = len( self.pdatafile.points )

        for i in range(niter):
            # use numarray arrays to store the data columns we are
            # interested in
            pgas = grabcolumn(self.pdatafile, i, 'pgas')
            pram = grabcolumn(self.pdatafile, i, 'pram')
            pmag = grabcolumn(self.pdatafile, i, 'P(mag)')
            vel = grabcolumn(self.pdatafile, i, 'windv')
            vel = abs(vel)*1.e5 # put in cgs units (and positive)
            hden = grabcolumn(self.odatafile, i, 'hden')
            hden = self.amass*10**(hden) # put in g/cm^3
            dr = grabcolumn(self.drdatafile, i, 'dr')
            depth = numarray.cumsum(dr)
            try: 
                # specific enthalpies
                hram = 0.5*pram/hden
                hgas = 2.5*pgas/hden
                hmag = 4.0*pmag/hden
                h = hram + hgas + hmag
                massflux = hden*vel
            except ValueError:
                # probably incompatible shapes due to a truncated file
                print 'Warning: incompatible array shapes'
                hram = numarray.zeros(len(dr),numarray.Float)
                hgas = numarray.zeros(len(dr),numarray.Float)
                hmag = numarray.zeros(len(dr),numarray.Float)
                h = numarray.zeros(len(dr),numarray.Float)
                massflux = numarray.ones(len(dr),numarray.Float)*hden[0]*vel[0]

            # calculate gradient of enthalpy by fitting straight line
            # to (2n+1)-point moving window
            dhdz = findgradient( depth, h, n )

            self.titles = [ "depth", "hram", "hgas", "hmag", "h", 
			    "rho*v", "dh/dz" ]

            # fill a table of rows for this iteration
            try: 
                iteration = [[depth[j], hram[j], hgas[j], hmag[j], 
			      h[j], massflux[j], dhdz[j]] 
			     for j in range(len(dr))]
            except IndexError:
                # probably one of the arrays is too short (incomplete file)
                iteration = []
            # now stuff into the points list
            self.points.append(iteration)
    
class heatcoolfile(dumbdatafile) :
    """Reads a heating or cooling file, which consist of a list of the
    important processes at each depth. Try to generatea a column of
    rates for each coolant with zeros inserted where the coolant is
    not signifiacant."""

    def __init__(self, file):
        "Read heating or cooling data from a file"
        self.points = []
        self.nzones = []
        dumbdatafile.__init__(self,file)
        for iteration in self.iterations:
            rows = iteration.split(self.row_separator)
            if not self.titles:
                # take the column headers from the first line of first
                # iteration after removing the `#'
                self.titles  = rows[0][1:].split(self.field_separator)
                # eliminate the last column (heat fracs or cool fracs)
                self.titles.remove(self.titles[-1])
                # number of fixed columns
                nfixed = len(self.titles)

            # each iteration is overwritten as a dictionary of form
            # {header0: data0, ...}
            iteration = {} 
            # guess at number of lines (probably overestimate)
            nrows = len(rows)
            # add dict entries for each title
            for title in self.titles:
                    iteration[title] = numarray.zeros(nrows,numarray.Float)

            irow = 0
            for row in rows:
                # First three fields are fixed
                # Subsequent fields are 'coolant\t value'
                try: 
                    fields = row.split(self.field_separator)
                    # first the fixed columns
                    for title in self.titles:
                        # set each fixed field
                        iteration[title][irow] = maybefloat(fields[0])
                        # and remove it from the field
                        fields.remove(fields[0])
                        # now the remaining fields (heat/cool agents)
                    while fields :
                        value, key = maybefloat(fields.pop()), fields.pop()
                        try:
                            # see if we already have a column for this
                            # coolant...
                            iteration[key][irow] = value
                        except KeyError:
                            # ...otherwise create one first
                            iteration[key] = numarray.zeros(nrows,numarray.Float)
                            iteration[key][irow] = value
                        # assume second column is the total heating
                        # and multiply each coolant by that
                        iteration[key][irow] *= iteration[self.titles[1]][irow]
                    irow += 1 # increment row counter
                except:
                    # Failed for some reason - probably an empty row
                    # Just ignore it
                    continue
            
            # now that this iteration is sorted, we can stuff it into points
            self.points.append( iteration )
            # and save the actual number of points
            self.nzones.append(irow)
            
    def writefile(self, prefix):
        i = 0
        for iteration in self.points:
            filename = prefix + "%4.4i" % (i+1)
            try:
                f = open(datadir + filename, "w")
            except:
                print "Error opening cloudy data file: %s" % filename
            keys = iteration.keys()
            sortkeys = []
            # rearrange to have the fixed columns first
            for key in self.titles:
                sortkeys.append(key)
                keys.remove(key)
            # and advection next (if present)
            for key in keys:
                if key.startswith("adve"):
                    sortkeys.append(key)
                    keys.remove(key)
            # and the rest in alphabetical order
            keys.sort()
            try: 
                sortkeys.extend(keys)
            except TypeError:
                print "No more keys for iteration %i" % (i)
            # write the header
            f.write("#" + "\t".join(sortkeys) + "\n")
            # write the data
            for jrow in range(self.nzones[i]):
                f.write("\t".join(
			["%5e" % iteration[key][jrow] 
			 for key in sortkeys]) + "\n")
            f.close()
            i += 1


                        




        
