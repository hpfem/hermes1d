from pylab import *
import numpy
import scipy

data = numpy.loadtxt("solution.gp_0")
U_re = data[:, 1]
x = data[:, 0]
data = numpy.loadtxt("solution.gp_1")
U_im = data[:, 1]

data = numpy.loadtxt("solution.gp_2")
x = data[:, 0]
I_re = data[:, 1]
data = numpy.loadtxt("solution.gp_3")
I_im = data[:, 1]
figure(1)
U = scipy.sqrt(U_re**2+U_im**2);
plot(x, U)
xlabel('l[m]');
ylabel('U[V]');
ylim(0,2)

I=scipy.sqrt(I_re**2+I_im**2);
figure(2)
plot(x, I);
xlabel('l[m]');
ylabel('I[A]');
ylim(0,0.035)
show();