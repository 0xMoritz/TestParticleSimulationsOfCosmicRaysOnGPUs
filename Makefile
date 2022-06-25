# Path to Cuda Installation
CUDA_ROOT = /usr/local/cuda-11.1

# target architecture
ARCH = sm_61

# Common Flags (-03 specifies level of optimization)
COMMON_FLAGS = -O3 -std=c++14

# Path to NVCC installation
NVCC = $(CUDA_ROOT)/bin/nvcc

# NVCC Flags
NVCC_FLAGS = $(COMMON_FLAGS) -arch $(ARCH) -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored


# Gnu C++ Compiler
#GCC = g++

# GCC Flags
#GCC_FLAGS = $(COMMON_FLAGS)

# All *.cpp files and all *.cu files
TARGETS = $(addsuffix .o, $(basename $(wildcard *.cpp))) $(addsuffix .o, $(basename $(wildcard *.cu)))

# Executable
BINARIES = CosmicRayTransport.bin


all : $(BINARIES)


%.o : %.cpp
	@echo "Building " $@
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

%.o : %.cu
	@echo "Building " $@
	$(NVCC) $(NVCC_FLAGS) -o $@ -c $<

CosmicRayTransport.bin : $(TARGETS)
	@echo "Building " $@
	$(NVCC) $(NVCC_FLAGS) $^ -o $@
	@echo "Finished Building."
	@echo ""

PHONY. : clean

clean:
	rm -f $(TARGETS)
	rm -f $(BINARIES)

