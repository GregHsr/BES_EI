# Nom du compilateur
FC = gfortran

# Options de compilation : optimisation, debug, etc.
OPT = -fdefault-real-8 -fcheck=all -g
 #-O3

# Nom de l'exécutable
EXE = cavite.exe

# Variables pour PETSc
PETSC_DIR = /home/greg/install/petsc      
PETSC_ARCH = arch-linux-c-debug 

# Inclure les fichiers d'en-tête PETSc
INCLUDE = -I/home/greg/install/petsc/include -I/home/greg/install/petsc/arch-linux-c-debug/include

# Lien avec les bibliothèques PETSc
LIBS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc

# Options d'édition de lien (ajout des bibliothèques PETSc)
LINKOPT = $(LIBS)

# Fichiers objets
OBJS =  \
       main.o \
       matgen_cavite.o \
       Solveur.o \
       ensight_cavite.o \
       advection.o

# Cible principale : génération de l'exécutable
exe : $(OBJS)
	$(FC) $(LINKOPT) $(OBJS) -o $(EXE)

# Règles pour compiler les fichiers sources
main.o : main.f90
	$(FC) -c $(OPT) $(INCLUDE) main.f90

matgen_cavite.o : matgen_cavite.f90
	$(FC) -c $(OPT) $(INCLUDE) matgen_cavite.f90

Solveur.o : Solveur.f90
	$(FC) -c $(OPT) $(INCLUDE) Solveur.f90

ensight_cavite.o : ensight_cavite.f90
	$(FC) -c $(OPT) $(INCLUDE) ensight_cavite.f90

advection.o : advection.f90
	$(FC) -c $(OPT) $(INCLUDE) advection.f90

# Cible pour nettoyer les fichiers objets et modules
clean :
	rm -f $(OBJS) $(EXE) *.mod

# Configuration : création des répertoires pour les fichiers objets et d'exécution
config :
	mkdir -p obj run
