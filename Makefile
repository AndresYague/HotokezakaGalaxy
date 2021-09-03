#
# Makefile
#
# snuppat
#

# Compiler.
FC = mpifort

# Flags.
 
# Optimization flags.
FFLAGS = -march=native -O2
MODFLAGS = -c -march=native -O2

# Profile flags.
ifeq ($(mode), p)
    FFLAGS = -march=native -O2 -g
    MODFLAGS = -c -march=native -O2 -g
endif

# Debug flags.
ifeq ($(mode), d)
    FFLAGS = -Og -ggdb3 -Wall -Wextra -Wconversion -pedantic -fbacktrace\
			 -ffpe-trap=invalid,zero,overflow -fbounds-check
    MODFLAGS = -c -Og -ggdb3 -Wall -Wextra -Wconversion -pedantic -fbacktrace\
			   -ffpe-trap=invalid,zero,overflow -fbounds-check
endif

# This Makefile
MAKEFILE = Makefile

# Executables.
EXE = HotokezakaCalc

# All:
ALL = $(EXE) tags

#--------------------------------------------

# All rule
all: $(ALL) $(MAKEFILE)

# Rules

# Ctags
tags: *.f90 $(MAKEFILE)
	ctags *.f90

# Executable compilation.
%: %.f90
	$(FC) $(FFLAGS) $(filter-out $(MAKEFILE), $^) -o $@

# Cleaning.
clean:
	rm -f $(filter-out tags, $(ALL)) $(MODS:=.mod) $(OBJECTS:=.mod)
