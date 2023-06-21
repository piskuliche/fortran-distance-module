# Compiler settings
FC = gfortran
FFLAGS = -O3 -std=f2008 

# Source files
SRC = main.f90 distance_module.f90 testing.f90

# Object files
OBJ = $(SRC:.f90=.o)

# Executable name
EXE = my_program

# Targets
all: $(EXE)

$(EXE): $(OBJ)
    $(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
    $(FC) $(FFLAGS) -c $< -o $@

clean:
    rm -f $(OBJ) $(EXE)