CC=g++
NVCC=nvcc
CFLAGS=-O3 -Wall
NVCCFLAGS=-O3 --use_fast_math -Wall

CC_LIBS=
NVCC_LIBS=

# CUDA library directory:
CUDA_LIB_DIR= -L$(CUDA_ROOT_DIR)/lib64
# CUDA include directory:
CUDA_INC_DIR= -I$(CUDA_ROOT_DIR)/include

INCLUDE_DIR=include
SRC_DIR=src
OBJ_DIR=bin



EXE=cif2ztb

OBJS = $(OBJ_DIR)/main.o $(OBJ_DIR)/energy.o $(OBJ_DIR)/cif2ztb_kernel.o $(OBJ_DIR)/unitcell.o $(OBJ_DIR)/symmetry.o

$(EXE) : $(OBJS)
	$(NVCC) $(NVCC_FLAGS) $(OBJS) -o $@ $(CUDA_INC_DIR) $(CUDA_LIB_DIR) -I $(INCLUDE_DIR)

# Compile main .cpp file to object files:
$(OBJ_DIR)/main.o :  $(SRC_DIR)/main.cpp 
	$(CC) $(CC_FLAGS) -c $< -o $@ -I $(INCLUDE_DIR)

# Compile CUDA source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cu
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(NVCC_LIBS) -I $(INCLUDE_DIR)

# Compile C++ source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CC) $(CC_FLAGS) -c $< -o $@ -I $(INCLUDE_DIR)


# Clean objects in object directory.
clean:
	$(RM) bin/* *.o $(EXE)
