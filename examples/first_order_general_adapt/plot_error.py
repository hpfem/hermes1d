from pylab import plot, show, legend
import numpy
data = numpy.loadtxt("error.gp")
x = data[:, 0]
y = data[:, 1]
plot(x, y, label="error")
legend()
show()
