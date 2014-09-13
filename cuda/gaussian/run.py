#!us/bin/python
#FILENAME:run.py

import os


os.system("nvcc track.cu -I. -I/usr/local/cuda/include  -L/usr/local/cuda/lib64 -arch=sm_21 `gsl-config --libs`   -O3 -w")
os.system("./a.out")