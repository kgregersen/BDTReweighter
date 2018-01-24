# Directories
OBJ=obj
SRC=src
INC=inc
BIN=bin


# Set ROOT variables
ROOTC = $(shell root-config --cflags)
ROOTLIB := $(shell root-config --libs)

# Set compiler flags
GCC = g++ -Wall -Wformat=0 -std=c++11
COPT = $(ROOTC) -I$(INC)


# Set linker flags
# Ubuntu's new default linker setting (--as-needed) exposes that many 
# libraries are linked incorrectly in ROOT. To fix this, we use the flag
# --no-as-needed for Ubuntu. Currently, this is not an issue for other 
# Linux distributions.
LD = g++
UNAME_OS := $(shell lsb_release -si)
ifeq ($(UNAME_OS),Ubuntu)
	LDFLAGS	= "-Wl,--no-as-needed" $(ROOTLIB) -L$(OBJ) # Ubuntu
else 
	LDFLAGS	= $(ROOTLIB) -L$(OBJ) # Other OS
endif


# List of sources and objects for main program
CXXSRC=$(shell find $(SRC) -name "*.cxx")
CXXOBJ=$(CXXSRC:$(SRC)/%.cxx=$(OBJ)/%.o)
CPPSRC=$(shell find $(SRC) -name "*.cpp")
CPPOBJ=$(CPPSRC:$(SRC)/%.cpp=$(OBJ)/%.o)
CPPEXE=$(CPPSRC:$(SRC)/%.cpp=$(BIN)/%)

# Set default target
all: $(CPPEXE)


# Generic rule for CXXOBJ
$(OBJ)/%.o: $(SRC)/%.cxx $(INC)/%.h 
	@echo " "
	@echo "------>>>>>> Compiling $<"
	$(GCC) $(COPT) -c $< -o $@

# Generic rule for CPPOBJ
$(OBJ)/%.o: $(SRC)/%.cpp
	@echo " "
	@echo "------>>>>>> Compiling $<"
	$(GCC) $(COPT) -c $< -o $@

# Link executables
$(BIN)/%: $(OBJ)/%.o $(CXXOBJ) 
	@echo " "
	@echo "------>>>>>> Linking $<"
	$(LD) $(LDFLAGS) $^ -o $@

# Clean
clean:
	@echo " "
	@echo "------>>>>>> Removing object files and executable"
	rm -rf $(BIN)/* $(OBJ)/*.o 

