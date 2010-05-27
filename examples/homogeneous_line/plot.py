import pylab
import numpy

N = 10
steps = 510
tu = numpy.empty([N, steps])
ti = numpy.empty([N, steps])
uu = numpy.empty([N, steps])
ii = numpy.empty([N, steps])

for i in range(1, N):
    data = numpy.loadtxt("solution_" + str(i-1) + ".gp_0")
    x = data[:, 0]
    y = data[:, 1]
    tu[i, :] = x
    uu[i, :] = y
    
    data = numpy.loadtxt("solution_" + str(i-1) + ".gp_1")
    x = data[:, 0]
    y = data[:, 1]
    ti[i, :] = x
    ii[i, :] = y

pylab.subplot(2,1,1)
pylab.plot(tu[:,1], uu[:,1], label="voltage")
pylab.grid(True);
pylab.title("voltage")

pylab.subplot(2,1,2)
pylab.plot(ti[:,0], ii[:,0], label="current")
pylab.grid(True);
pylab.title("current")
pylab.show()
