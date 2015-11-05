#!/usr/bin/env python2.7
from mpi4py import MPI
import random
import numpy as np

start_time = MPI.Wtime()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
mpisize = comm.Get_size()
nsamples = int(12e7/mpisize)
inside = 0

np.random.seed(rank)
xy=np.random.random((nsamples,2))
mypi=4.0*np.sum(np.sum(xy**2,1)<1)/nsamples

pi = comm.reduce(mypi, op=MPI.SUM, root=0)

end_time = MPI.Wtime()

if rank==0:
    print (1.0 / mpisize)*pi
    print "Run time for %s was %s" %(nsamples,end_time-start_time)
