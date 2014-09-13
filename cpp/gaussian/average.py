#!/usr/bin/python

import linecache
import math
import os



f1=open("abc.txt","r")
x=f1.readlines()

count=0
ex=0
delta=0
sdelta=0
for i in range(1700,2000):

	ex=ex+float(x[i][4:19])
#	delta=delta+float(x[i][34:48])
#	sdelta=sdelta+ (float(x[i][34:48])-8.12E-5) * ( float(x[i][34:48]) -8.12E-5)
	count=count+1
	

print ex/count
