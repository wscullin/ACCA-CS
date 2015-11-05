#!/usr/bin/env python2.7

import numpy as np
import time


start_time = time.time()

nsamples = int(12e7)
inside = 0

np.random.seed(np.int(start_time))
xy=np.random.random((nsamples,2))
pi=4.0*np.sum(np.sum(xy**2,1)<1)/nsamples

end_time = time.time()


print pi
print "Run time for %s was %s" %(nsamples,end_time-start_time)
