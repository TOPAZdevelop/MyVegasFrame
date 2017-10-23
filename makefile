Here = $(PWD)
ConfigFile = $(Here)/config.cfg
ModuleDir = $(Here)/modules
ObjectDir = $(Here)/objects
VegasDir = $(Here)/Vegas
PSDir = $(Here)/PhaseSpace
QCDLoop = $(Here)/QCDLoop-1.9

# For ifort support on CERN's lxplus, add the following two lines to your .bashrc
# source /afs/cern.ch/sw/IntelSoftware/linux/x86_64/xe2016/compilers_and_libraries_2016.1.150/linux/bin/ifortvars.sh intel64
# source /afs/cern.ch/sw/IntelSoftware/linux/setup.sh


# ifort optimization, Yes/No
Opt = Yes

# MPI features, Yes/No. run with: mpiexec -n 4 ./NumInt_MPI ...
useMPI = No


ifeq ($(useMPI),Yes)
    Exec = ./NumInt_MPI
    F95compiler = mpif90-ifort -lpthread -D_UseMPIVegas=1
    ccomp = mpicc -lpthread  -lm 
else
    Exec = ./NumInt
    F95compiler = ifort -D_UseMPIVegas=0
    ccomp = gcc -O2
endif




ifeq ($(Opt),Yes)
   IfortOpts   = -O2 -fpp -I$(Here)/colors -I$(VegasDir) -module $(ModuleDir) 
else
   IfortOpts   = -O0 -fpp -implicitnone -check bounds -check pointer -warn interfaces -ftrapuv  -debug extended -g -traceback -fpe0 -check uninit -I$(Here)/colors -I$(VegasDir) -module $(ModuleDir) 
endif
fcomp = $(F95compiler) $(IfortOpts) @$(ConfigFile)






# makeDep = $(ConfigFile) makefile


# fastjet stuff
#FASTJET_CONFIG=/home/schulze/usr/local/bin/fastjet-config
# CXXFLAGS += $(shell $(FASTJET_CONFIG) --cxxflags)
# FJLIBS += $(shell $(FASTJET_CONFIG) --libs --plugins )


         
RockyObj = $(ObjectDir)/genps.o \
           $(ObjectDir)/boost.o
           
ifeq ($(useMPI),Yes)
   VegasObj = $(VegasDir)/pvegas_mpi.o
else

   VegasObj = $(VegasDir)/nvegas.o
endif



IntegralObj = $(QCDLoop)/ql/libqcdloop.a\
              $(QCDLoop)/ff/libff.a




# ------------------------------------------------------------


# the order of these object files corresponds to their mutual dependency
allObjects =   				$(ObjectDir)/mod_Misc.o \
					$(ObjectDir)/mod_Parameters.o \
					$(ObjectDir)/mod_Amplitudes.o \
					$(ObjectDir)/mod_Kinematics.o \
					$(ObjectDir)/mod_Integrand.o \
					$(ObjectDir)/main.o





all:  $(VegasObj) $(RockyObj) $(allObjects)
	@echo " linking"
	@echo " executable file is " $(Exec)
	@echo " "
	$(fcomp) -o $(Exec) $(allObjects) $(RockyObj) $(VegasObj) 

clean:
	rm -f ./modules/*.mod
	rm -f ./objects/*.o
	rm -f ./Vegas/*.o



./summer: summer.f90
	@echo " compiling" $<
	$(fcomp) $< -o summer


# ------------------------------------------------------------





$(ObjectDir)/mod_Misc.o: mod_Misc.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Parameters.o: mod_Parameters.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Amplitudes.o: mod_Amplitudes.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


	
$(ObjectDir)/main.o: main.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@


$(ObjectDir)/mod_Integrand.o: mod_Integrand.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@




$(ObjectDir)/mod_Kinematics.o: mod_Kinematics.f90 $(makeDep)
	@echo " compiling" $<
	$(fcomp) -c $< -o $@

$(VegasDir)/pvegas_mpi.o: $(VegasDir)/pvegas_mpi.c $(makeDep)
	@echo " compiling" $<
	$(ccomp) -D_WriteTmpHisto=1 -c $(VegasDir)/pvegas_mpi.c -o $@ 

$(VegasDir)/nvegas.o: $(VegasDir)/nvegas.c $(makeDep)
	@echo " compiling" $<
	$(ccomp) -D_WriteTmpHisto=1 -c $(VegasDir)/nvegas.c -o $@

	
$(ObjectDir)/genps.o: $(PSDir)/genps.c $(makeDep)
	@echo " compiling" $<
	$(ccomp) -c $< -o $@

	
$(ObjectDir)/boost.o: $(PSDir)/boost.c $(makeDep)
	@echo " compiling" $<
	$(ccomp) -c $< -o $@

	

# supresses command calls
.SILENT:
