#!/usr/bin/env python2.7

import random 
import time


start_time = time.time()

nsamples = int(12e7)
inside = 0

random.seed(start_time)
pi=4.0*sum([random.random()**2+random.random()**2<1 for i in range(nsamples)])/nsamples

end_time = time.time()


print pi
print "Run time for %s was %s" %(nsamples,end_time-start_time)
