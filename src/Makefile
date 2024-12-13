# Makefile 
# Nom du compilateur
FC = gfortran

# Options de compilation: optimisation, debug etc...
OPT = -fdefault-real-8  -C -g
# nom de l'executable
EXE = cavite.exe
# Options de l'edition de lien..
LINKOPT =  

# Defining the objects (OBJS) variables

OBJS =  \
       main.o \
       matgen_cavite.o \
       Solveur.o \
       ensight_cavite.o \
	   advection.o

# Linking object files
exe :   $(OBJS)
	$(FC) $(LINKOPT) $(MODS) $(OBJS)  -o $(EXE) 

main.o : main.f90
	$(FC) -c $(OPT)  main.f90

matgen_cavite.o : matgen_cavite.f90
	$(FC) -c $(OPT) matgen_cavite.f90

Solveur.o : Solveur.f90
	$(FC) -c $(OPT) Solveur.f90

ensight_cavite.o : ensight_cavite.f90
	$(FC) -c $(OPT) ensight_cavite.f90

advection.o : advection.f90
	$(FC) -c $(OPT) advection.f90

# Removing object files
clean :
	/bin/rm -f $(OBJS) $(EXE)  *.mod

config :
	if [ ! -d obj ] ; then mkdir obj ; fi
	if [ ! -d run ] ; then mkdir run ; fi
