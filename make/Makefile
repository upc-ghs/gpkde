# makefile created on 2017-12-20 13:16:47.325722
# by pymake (version 1.1.0)
# using the gfortran fortran and gcc c/c++ compilers.

# Define the directories for the object and module files,
# the executable, and the executable name and path.
OBJDIR = ./objtemp
BINDIR = .
PROGRAM = gpkde

SOURCEDIR1=../src

VPATH = \
${SOURCEDIR1}


.SUFFIXES: .c .cpp .f .f90 .F90 .fpp .o

# Define the Fortran compile flags
FC = gfortran
#FFLAGS = -O2 -g -pg -fbacktrace -fopenmp -Wall -Wextra -fbounds-check -fcheck=bounds -fcheck=all
#FFLAGS =  -g -pg -fbacktrace -fopenmp
FFLAGS = -O2 -g -pg -fbacktrace -fopenmp
#FFLAGS = -g -pg -fbacktrace -fopenmp -Wextra -fcheck=bounds -fcheck=all
#FFLAGS = -O2 -fbacktrace -fopenmp

# Define the C compile flags
CC = gcc
CFLAGS = -O3 -D_UF 

# Define the libraries
SYSLIBS = 

OBJECTS = \
$(OBJDIR)/Histogram.o \
$(OBJDIR)/KernelMultiGaussian.o \
$(OBJDIR)/GridCell.o \
$(OBJDIR)/GridProjectedKDE.o \
$(OBJDIR)/testkerneldbreal.o
#$(OBJDIR)/testgpkdeomp.o
#$(OBJDIR)/testgpkde.o
#$(OBJDIR)/testunfoldmatrix.o

# Define task functions

# Create the bin directory and compile and link the program
all: makebin | $(PROGRAM)

# Make the bin directory for the executable
makebin :
	mkdir -p $(BINDIR)

# Define the objects that make up the program
$(PROGRAM) : $(OBJECTS)
	-$(FC) $(FFLAGS) -o $@ $(OBJECTS) $(SYSLIBS) -I$(OBJDIR) -J$(OBJDIR) 

$(OBJDIR)/%.o : %.f
	@mkdir -p $(@D)
	$(FC) $(FFLAGS) -c $< -o $@ -I$(OBJDIR) -J$(OBJDIR) 

$(OBJDIR)/%.o : %.f90
	@mkdir -p $(@D)
	$(FC) $(FFLAGS) -c $< -o $@ -I$(OBJDIR) -J$(OBJDIR) 

$(OBJDIR)/%.o : %.F90
	@mkdir -p $(@D)
	$(FC) $(FFLAGS) -c $< -o $@ -I$(OBJDIR) -J$(OBJDIR) 

$(OBJDIR)/%.o : %.fpp
	@mkdir -p $(@D)
	$(FC) $(FFLAGS) -c $< -o $@ -I$(OBJDIR) -J$(OBJDIR) 

$(OBJDIR)/%.o : %.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o : %.cpp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean the object and module files and the executable
.PHONY : clean
clean : 
	-rm -rf $(OBJDIR)
	-rm -rf $(PROGRAM)

# Clean the object and module files
.PHONY : cleanobj
cleanobj : 
	-rm -rf $(OBJDIR)

