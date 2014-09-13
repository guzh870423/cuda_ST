
CUDA_INSTALL_PATH ?= /usr/local/cuda
NVCC       := $(CUDA_INSTALL_PATH)/bin/nvcc 
INCLUDES  += -I. -I$(CUDA_INSTALL_PATH)/include 
LINK := -L$(CUDA_INSTALL_PATH)/lib64
file :=track.cpp
ARCH := -arch=sm_21
GSL := `gsl-config --libs`
BOOST := -I /home/guzh/boost_1_55_0/


a.out   :   $(file)
	g++ $(file) $(GSL) $(BOOST) -lm -O3 -w

