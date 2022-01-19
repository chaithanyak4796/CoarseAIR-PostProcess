import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import sys

# plt.close('all')

au_s = 2.4188843265E-17;
Bo_angs = 0.529177

if(len(sys.argv) == 3):
    iP    = int(sys.argv[1])
    iTraj = int(sys.argv[2])
else:
    iP = 1
    iTraj = 36

Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Test_Current/T_8000_8000_0_0/Bins_1951_0/Node_1/Proc_"+str(iP)

fname = Dir + "/PaQEvo-" + str(iTraj) + ".out"

data = np.loadtxt(fname,skiprows=1)
t = data[:,0]
t = t*au_s*1E12

Q = np.zeros([len(t),9])
P = np.zeros([len(t),9])

P[:,0:6] = data[:,2:8] /(au_s*1E12)
Q[:,0:6] = data[:,8:14]

Q[:,6] = - Q[:,0] - Q[:,3]
Q[:,7] = - Q[:,1] - Q[:,4]
Q[:,8] = - Q[:,2] - Q[:,5]

P[:,6] = - P[:,0] - P[:,3]
P[:,7] = - P[:,1] - P[:,4]
P[:,8] = - P[:,2] - P[:,5]

R = np.zeros([len(t),4])
rho = np.zeros((len(t),))
E = np.zeros((len(t),2))


R[:,0] = ( (Q[:,0]-Q[:,3])**2 + (Q[:,1]-Q[:,4])**2 + (Q[:,2]-Q[:,5])**2 )**0.5
R[:,1] = ( (Q[:,0]-Q[:,6])**2 + (Q[:,1]-Q[:,7])**2 + (Q[: ,2]-Q[:,8])**2 )**0.5
R[:,2] = ( (Q[:,3]-Q[:,6])**2 + (Q[:,4]-Q[:,7])**2 + (Q[:,5]-Q[:,8])**2 )**0.5
R[:,3] = ( (Q[:,6]-0.5*(Q[:,0]+Q[:,3]))**2 + (Q[:,7]-0.5*(Q[:,1]+Q[:,4]))**2 + (Q[:,8]-0.5*(Q[:,2]+Q[:,5]))**2 )**0.5
# rho = (R[:,0]**2 + R[:,1]**2 + R[:,2]**2)**0.5
rho = (R[:,0]**2 + R[:,3]**2 )**0.5

plt.figure(1)
lt = '-'
plt.plot(t,R[:,0],lt,label=r'R$_{AB}$')
plt.plot(t,R[:,1],lt,label=r'R$_{AC}$')
plt.plot(t,R[:,2],lt,label=r'R$_{BC}$')
# plt.plot(t,R[:,3],lt,label=r'R$_{AB-C}$')
# plt.plot(t,rho,lt,label=r'$\rho$')
plt.xlabel('time [ps]')
plt.ylabel(r'Inter-nuclear distance [$\AA$]',fontweight='bold')
plt.legend()
#plt.show(block=False)
plt.show()
