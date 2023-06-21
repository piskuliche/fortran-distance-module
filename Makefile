# Compiler settings
FC = gfortran
FFLAGS = -O3 -std=f2008 -fbounds-check

# Source files
SRC_DIR = src
MOD = $(wildcard $(SRC_DIR)/*module.f90)
SRC = $(wildcard $(SRC_DIR)/*.f90)

# Object files
OBJ_DIR = obj
MOD_OBJ = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(MOD))
MAIN_OBJ = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(filter-out $(MOD),$(SRC)))

# Executable name
EXE = my_program

# Targets
all: $(EXE)

$(EXE): $(MOD_OBJ) $(MAIN_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm $(OBJ_DIR)/*.o $(EXE) *.mod

.PHONY: clean