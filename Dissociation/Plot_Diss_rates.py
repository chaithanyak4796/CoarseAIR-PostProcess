import numpy as np
import matplotlib.pyplot as plt

system = "N3"

if(system == "O3"):
    data = np.array([[8000.0,5.73619773E-12],
                     [6000.0,5.24267615E-13]],float)
    data_lit = np.loadtxt('O2_Diss.dat')
elif(system == "N3"):
    data = np.array([[20000.0,1.73893648E-11],
                     [18000.0,9.76632267E-12],
                     [15000.0,3.11228137E-12],
                     [13000.0,1.03172359E-12],
                     [10000.0,8.07347389E-14]],float)

    #data_lit = np.array([[8000,4.6E-15],[10000,7.1E-14],[13000,9.13E-13],[20000,1.633E-11]],float)
    data_lit = np.array([[10000,1.310E-13],[12500,1.081E-12],[15000,4.337E-12],[20000,2.280E-11]],float)
    data_lit[:,0] = 10000/data_lit[:,0]

plt.semilogy(data_lit[:,0],data_lit[:,1],label='Literature')
plt.semilogy(10000/data[:,0],data[:,1],'o',label='My Calcs')
plt.xlabel(r'$10000/T [K^{-1}]$')
plt.ylabel(r'$k_{diss} [cm^3/s]$')
plt.legend()
plt.show()
