include arch.mk

EXECS=cubeextras cubet enro2fits cubeemiss makerotmap makerotbmaps globstats

all: $(EXECS)

clean: 
	rm -fr $(EXECS) *.o *.mod

SOURCES=*.{f90,py,sh} Makefile arch-*.mk

THISDIR=Fabio/PostProcess/src
ship2pinga:
	rsync -avzuP $(SOURCES) \
              pinga:/fs/pinga/other0/will/$(THISDIR)
ship2magneton:
	rsync -avzuLP $(SOURCES) \
              magneton:/fs/magneton/other0/will/$(THISDIR)
ship2manouapa:
	rsync -avzuLP $(SOURCES) \
              manouapa:/fs/manouapa/other0/will/$(THISDIR)
ship2robbie:
	rsync -avzuLP $(SOURCES) \
              robbie:/fs/robbie/other0/will/$(THISDIR)
ship2kanbalam:			# for kanbalam, copy target of links since dir layout is different 
	rsync -avzuLP $(SOURCES) \
              kanbalam:/global/home/whenney_g/whenney/$(THISDIR)
ship2pinga_remote:		# This requires the relevant tunnel to be operating
	rsync -e "ssh -p 2222" -avzuP $(SOURCES) \
              localhost:/fs/pinga/other0/will/$(THISDIR)

cubeextras.o:  wfitsutils.o
cubeemiss.o:  wfitsutils.o emissmod.o em2levmod.o
emissmod.o: em2levmod.o
rotatetest.o: cuberotate.o wfitsutils.o
makerotmap.o: cuberotate.o wfitsutils.o emissmod.o
makerotbmaps.o: cuberotate.o wfitsutils.o molfrac.o
cubet.o:  wfitsutils.o
enro2fits.o:  wfitsutils.o

globcollate.o: allstats_vars.o
globstats.o:  wfitsutils.o modstats.o allstats_vars.o
turbstats.o:  wfitsutils.o turbstats_vars.o
cubestats.o:  wfitsutils.o
cubevstats.o:  wfitsutils.o
cubetstats.o:  wfitsutils.o
cubedenstats.o:  wfitsutils.o
omp_test.o omp_test2.0 omp_test3.o omp_test4.o:  wfitsutils.o

phabrebin.o: wfitsutils.o

phabrebin: phabrebin.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

cubestats:    cubestats.o    wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
cubevstats:   cubevstats.o   wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
cubetstats:   cubetstats.o   wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
cubedenstats: cubedenstats.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)


turbstats: turbstats.o wfitsutils.o turbstats_vars.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
globstats: globstats.o wfitsutils.o modstats.o allstats_vars.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
globcollate: globcollate.o allstats_vars.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

cubeemiss: cubeemiss.o wfitsutils.o emissmod.o em2levmod.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
cubeextras: cubeextras.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
rotatetest: rotatetest.o cuberotate.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
makerotmap: makerotmap.o cuberotate.o wfitsutils.o emissmod.o em2levmod.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
makerotbmaps: makerotbmaps.o cuberotate.o wfitsutils.o molfrac.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

cubet: cubet.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

enro2fits: enro2fits.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

omp_test: omp_test.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
omp_test2: omp_test2.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
omp_test3: omp_test3.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
omp_test4: omp_test4.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
