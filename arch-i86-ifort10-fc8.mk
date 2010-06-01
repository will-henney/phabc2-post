MAKEDIR=$(HOME)/MAKE
MACHFILE=$(MAKEDIR)/linux-ifort.mk

include $(MACHFILE)

F90FLAGS=-xW -O3 -vec_report -u -ip -ipo -fpe0
# F90FLAGS=$(F90FLAGS_DEBUG)	
F77FLAGS=$(F90FLAGS) 

LFITS=-L$(HOME)/local/lib -lcfitsio
#LFITS=/usr/lib64/libcfitsio.so.0

FPP=-fpp

