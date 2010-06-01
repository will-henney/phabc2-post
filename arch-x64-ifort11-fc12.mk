MAKEDIR=$(HOME)/MAKE
MACHFILE=$(MAKEDIR)/linux-ifort-64.mk

include $(MACHFILE)

# version of -fast that doesn't use -static
FAST=-O3 -no-prec-div -xHost -ipo
UNSAFE=-fp-model fast=2 -fp-speculation=fast
PROFILE=-pg
OPENMP=-openmp -openmp-report2
F90FLAGS=-fpe0 $(FAST) $(UNSAFE) -vec_report -u $(OPENMP)
# F90FLAGS=-g -CB -debug extended -O0 -fpe0 -u
# F77FLAGS=$(F90FLAGS) -I $(INCLUDES) -w=unused

#LFITS=-Bdynamic -L/usr/lib64 -lcfitsio
LFITS=-Bdynamic /usr/lib64/libcfitsio.so.0

FPP=-fpp

