FC = ifort

##########################
clm.exe: clm_main.o shr_const_mod.o module_biogeophysics.o clmdata.o \
		module_initialize_clm.o
	$(FC) -o clm.exe clm_main.o module_biogeophysics.o clmdata.o \
		module_initialize_clm.o

clm_main.o: clm_main.f90 shr_const_mod.o module_biogeophysics.o clmdata.o \
		module_initialize_clm.o
	$(FC) -c clm_main.f90

module_biogeophysics.o: module_biogeophysics.f90 shr_const_mod.o clmdata.o
	$(FC) -c -cpp module_biogeophysics.f90

shr_const_mod.o: module_shr_const.f90
	$(FC) -c module_shr_const.f90

clmdata.o: clmdata.f90 shr_const_mod.o
	$(FC) -c -cpp clmdata.f90

module_initialize_clm.o: module_initialize_clm.f90 clmdata.o module_biogeophysics.o
	$(FC) -c  module_initialize_clm.f90

#module_canopy.o: module_canopy.f90 shr_const_mod.o clmdata.o
#	$(FC) -c module_canopy.f90

#####################

clean: 
	rm -f clm.exe *.o *.mod

