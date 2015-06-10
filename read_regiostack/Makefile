# $Id: Makefile,v 1.2 2005/06/19 20:29:35 alessia Exp $

OPT= -Wall -O2 -fbounds-check 
CC = gcc
FC = gfortran
CFLAGS = ${OPT}
FFLAGS = ${OPT}


BINDIR= ${HOME}/bin

SUFFIXES = .f90 .o
.SUFFIXES: $(SUFFIXES)

.f90.o:
	${FC} $(FFLAGS) -c -o $@ $<
	

depth_cut_iso: constants.h depth_cut_iso.f90 mantle_model_regiostack_iso.o  
	${FC} $(FFLAGS) -o depth_cut_iso depth_cut_iso.f90 mantle_model_regiostack_iso.o 

depth_cut_aniso: constants.h depth_cut_aniso.f90 mantle_model_regiostack_aniso.o  
	${FC} $(FFLAGS) -o depth_cut_aniso depth_cut_aniso.f90 mantle_model_regiostack_aniso.o 


gc_cut: constants.h gc_cut.f90 mantle_model_regiostack_iso.o 
	${F90} $(FLAGS_CHECK) -o gc_cut gc_cut.f90 mantle_model_regiostack_iso.o 

gc_cut_abs: constants.h gc_cut_abs.f90 mantle_model_regiostack_iso.o  
	${F90} $(FLAGS_CHECK) -o gc_cut_abs gc_cut_abs.f90 mantle_model_regiostack_iso.o 


get_moho: constants.h get_moho.f90 crustal_model_3smac.o rthetaphi_xyz.o reduce.o
	${F90} $(FLAGS_CHECK) -o get_moho get_moho.f90 crustal_model_3smac.o rthetaphi_xyz.o reduce.o

crust_cut: constants.h crust_cut.f90 crustal_model_3smac.o rthetaphi_xyz.o reduce.o
	${F90} $(FLAGS_CHECK) -o crust_cut crust_cut.f90 crustal_model_3smac.o rthetaphi_xyz.o reduce.o

clean:
	rm *.o
