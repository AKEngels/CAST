OBJS=configuration.o coords.o coords_io.o energy.o energy_int_aco.o energy_int_aco_pot.o energy_int_amoeba.o energy_int_amoeba_pot.o energy_int_mopac.o energy_int_terachem.o main.o matop.o md.o mpi_cast.o neb.o optimization.o optimization_grid.o optimization_mc.o optimization_ts.o Path_perp.o pathopt.o pme.o pme_fep.o reaccoord.o scon.o startopt.o startopt_ringsearch.o startopt_solvadd.o tinker_parameters.o tinker_refine.o

OMPBJS=configuration.om coords.om coords_io.om energy.om energy_int_aco.om energy_int_aco_pot.om energy_int_amoeba.om energy_int_amoeba_pot.om energy_int_mopac.om energy_int_terachem.om main.om matop.om md.om mpi_cast.om neb.om optimization.om optimization_grid.om optimization_mc.om optimization_ts.om Path_perp.om pathopt.om pme.om pme_fep.om reaccoord.om scon.om startopt.om startopt_ringsearch.om startopt_solvadd.om tinker_parameters.om tinker_refine.om

CC=g++
C=gcc
WARN=-Wextra -Wall -pedantic
DEB=
MPI=
CCSTD=-std=c++0x
CSTD=-std=c11
OPT=-O3
OPTCL=-flto
OMP=-fopenmp
LTO=-flto
STATIC=-static
EXE=C3.xgcc
PEXE=C3_omp.xgcc

VPATH=../src
BUILD=./
SOURCE=../src/


FLAGS=$(WARN) $(DEB) $(OPT) $(MPI) $(OPTCL)
CCFLAGS=-c $(CCSTD) $(FLAGS)
CFLAGS=-c $(CSTD) $(FLAGS)
LFLAGS=$(STD) $(OPTCL) $(STATIC)

CAST: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(EXE)

COMP: $(OMPBJS)
	$(CC) $(LFLAGS) $(OMP) $(OMPBJS) -o $(PEXE)

$(BUILD)%.o : %.cc
	$(CC) $(CCFLAGS) $<

$(BUILD)%.om : %.cc
	$(CC) $(OMP) $(CCFLAGS) $< -o $(<F:.cc=.om)

%.o : %.c
	$(C) $(CFLAGS) $<

%.om : %.c
	$(C) $(OMP) $(CFLAGS) $< -o $(<:.c=.om)

clean :
	rm $(OBJS) $(OMPBJS)

