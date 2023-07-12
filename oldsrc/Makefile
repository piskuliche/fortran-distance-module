# Compiler settings
FC = mpif90
FFLAGS = -O3 -std=f2008 -fopenmp -lmpi

# Source files
SRC_DIR = src
MOD = $(wildcard $(SRC_DIR)/*module.f90)
SRC = $(wildcard $(SRC_DIR)/*.f90)
MPI_MAIN_SRC = $(SRC_DIR)/mpi_main.f90

# Object files
OBJ_DIR = obj
MOD_SORTED = $(sort $(MOD))
MPI_MOD = $(filter %mpi_module.f90,$(MOD_SORTED))
NON_MPI_MOD = $(filter-out %mpi_module.f90,$(MOD_SORTED))
MOD_OBJ = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(NON_MPI_MOD)) $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(MPI_MOD))
MAIN_OBJ = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(filter-out $(MOD) $(MPI_MAIN_SRC),$(SRC)))
MPI_MAIN_OBJ = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(MPI_MAIN_SRC))

# Executable names
MAIN_EXE = my_program
MPI_MAIN_EXE = mpi_main

# Library names
LIB_DIR = lib
DISTANCE_LIB = $(LIB_DIR)/libdistance.so

# Targets
all: $(MAIN_EXE) $(MPI_MAIN_EXE) $(DISTANCE_LIB)

$(MAIN_EXE): $(MOD_OBJ) $(MAIN_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

$(MPI_MAIN_EXE): $(MOD_OBJ) $(MPI_MAIN_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(DISTANCE_LIB): $(OBJ_DIR)/distance_module.o
	$(FC) -shared -o $@ $<

clean:
	rm $(OBJ_DIR)/*.o $(MAIN_EXE) $(MPI_MAIN_EXE) $(DISTANCE_LIB) *.mod *.dat

.PHONY: clean