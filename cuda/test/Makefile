
CUDA_INSTALL_PATH ?= /usr/local/cuda
NVCC       := $(CUDA_INSTALL_PATH)/bin/nvcc 
INCLUDES  += -I. -I$(CUDA_INSTALL_PATH)/include 
LINK := -L$(CUDA_INSTALL_PATH)/lib64
file :=test.cu
ARCH := -arch=sm_21
GSL := `gsl-config --libs`

a.out : $(file)
	$(NVCC) $(file) $(INCLUDES) $(LINK) $(ARCH) $(GSL) -O3
