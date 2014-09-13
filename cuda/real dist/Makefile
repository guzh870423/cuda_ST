
CUDA_INSTALL_PATH ?= /usr/local/cuda
NVCC       := $(CUDA_INSTALL_PATH)/bin/nvcc 
INCLUDES  += -I. -I$(CUDA_INSTALL_PATH)/include 
LINK := -L$(CUDA_INSTALL_PATH)/lib64
file :=track.cu
ARCH := -arch=sm_21
GSL := `gsl-config --libs`
BOOST := -I /home/guzh/boost_1_55_0/
a.out :  $(file)
	$(NVCC) $(file) $(INCLUDES) $(LINK) $(ARCH) $(GSL)   -O3 -w
#a.out :   track.o  SynRad.o
#	$(NVCC) track.o SynRad.o  $(GSL)
#track.o   :   $(file)
#	$(NVCC) $(file) -c $(INCLUDES) $(LINK) $(ARCH) $(GSL)   -O3 -w
#SynRad.o  :  SynRad.cpp
#	g++ SynRad.cpp -c  $(GSL) -I /home/guzh/boost_1_55_0 -O3
