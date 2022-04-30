COMPILER=mpiifort
OPTIMIZATION=-Ofast
#VMC: io_handler.o bec_vmc.o VMC.f90
#	mpiifort -Ofast io_handler.o bec_vmc.o VMC.f90 -o VMC.x
#io_handler.o: io_handler.f90
#	mpiifort -Ofast -c io_handler.f90
#trial: bec_dmc.o trial.f90
#	$(COMPILER) -Ofast bec_dmc.o trial.f90 -o a.out
DMC: bec_dmc.o bec_vmc.o DMC.f90 
	$(COMPILER) -check all bec_dmc.o bec_vmc.o DMC.f90 -o DMC.x
bec_vmc.o: bec_vmc.f90 
	$(COMPILER) $(OPTIMIZATION) -c check all bec_vmc.f90	
bec_dmc.o: bec_dmc.f90 
	$(COMPILER) $(OPTIMIZATION) -c check all bec_dmc.f90
