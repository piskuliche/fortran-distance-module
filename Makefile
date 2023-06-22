# Compiler settings
FC = mpif90
FFLAGS = -O3 -std=f2008 -fopenmp -lmpi -fbounds-check

# Source files
SRC_DIR = src
MOD = $(wildcard $(SRC_DIR)/*module.f90)
SRC = $(wildcard $(SRC_DIR)/*.f90)
MPI_MAIN_SRC = $(SRC_DIR)/mpi_main.f90

# Object files
OBJ_DIR = obj
MOD_OBJ = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(MOD))
MAIN_OBJ = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(filter-out $(MOD) $(MPI_MAIN_SRC),$(SRC)))
MPI_MAIN_OBJ = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(MPI_MAIN_SRC))

# Executable names
MAIN_EXE = my_program
MPI_MAIN_EXE = mpi_main

# Targets
all: $(MAIN_EXE) $(MPI_MAIN_EXE)

$(MAIN_EXE): $(MOD_OBJ) $(MAIN_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

$(MPI_MAIN_EXE): $(MOD_OBJ) $(MPI_MAIN_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm $(OBJ_DIR)/*.o $(MAIN_EXE) $(MPI_MAIN_EXE) *.mod *.dat

.PHONY: clean