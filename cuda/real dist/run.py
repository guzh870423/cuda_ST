#!us/bin/python
#FILENAME:run.py

import os

#os.system("g++ SynRad.cpp `gsl-config --libs` -I'/home/guzh/boost_1_55_0/'")
#os.system("./a.out")
os.system("nvcc track.cu -I. -I/usr/local/cuda/include  -L/usr/local/cuda/lib64 -arch=sm_21 `gsl-config --libs`   -O3 ")
os.system("./a.out")