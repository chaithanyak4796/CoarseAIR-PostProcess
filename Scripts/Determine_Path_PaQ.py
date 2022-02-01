#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:30:15 2020

@author: ckondur
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import sys 

plt.close('all')

au_s = 2.4188843265E-17;
Bo_angs = 0.529177

if(len(sys.argv) == 3):
    iP    = int(sys.argv[1])
    iTraj = int(sys.argv[2])
else:
    iP = 4
    iTraj = 48

#Dir = "/Users/ckondur/Desktop/mount_sshfs/QCT_Output/Testing/O2/T_100_200_0_/Bins_10_0/Node_1/Proc_"+str(iP)

# Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Results/N3/Recombination/Rates/Pathway/Temp_18000K/T_18000_18000_0_11/Bins_10_0/Node_1/Proc_"+str(iP)
# Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/testing/CoarseAIR/N2O/Pathway/IncStpSz_yes/T_100_200_60_0/Bins_10_0/Node_1/Proc_" +str(iP)
Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Test_Current/T_100_200_60_0/Bins_10_0/Node_1/Proc_" +str(iP)

fname = Dir + "/PaQEvo-" + str(iTraj) + ".out"

too_far = False
invert_time = False

data = np.loadtxt(fname,skiprows=1)
if (len(np.array(np.shape(data))) == 1):
    data = data.reshape(1,-1) 
    too_far = True


t = data[:,0]
if(invert_time):
    t = t[-1] - t

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

t_unq = len(np.unique(t))
P_temp = [P[0]]
Q_temp = [Q[0]]
t_temp = [t[0]]
t_rep = []
idx_rep = []
for i in range(1,len(t)):
    if (t[i] > t_temp[-1]):
        t_temp.append(t[i])
        P_temp.append(P[i])
        Q_temp.append(Q[i])
    else:
        t_rep.append(t[i]/au_s*1E-12)
        idx_rep.append(i)
        
P = np.array(P_temp)
Q = np.array(Q_temp)
t = np.array(t_temp)

print("Number of repitions = ",len(data)-len(t))

R = np.zeros([len(t),4])
rho = np.zeros((len(t),))
E = np.zeros((len(t),2))


R[:,0] = ( (Q[:,0]-Q[:,3])**2 + (Q[:,1]-Q[:,4])**2 + (Q[:,2]-Q[:,5])**2 )**0.5
R[:,1] = ( (Q[:,0]-Q[:,6])**2 + (Q[:,1]-Q[:,7])**2 + (Q[: ,2]-Q[:,8])**2 )**0.5
R[:,2] = ( (Q[:,3]-Q[:,6])**2 + (Q[:,4]-Q[:,7])**2 + (Q[:,5]-Q[:,8])**2 )**0.5
R[:,3] = ( (Q[:,6]-0.5*(Q[:,0]+Q[:,3]))**2 + (Q[:,7]-0.5*(Q[:,1]+Q[:,4]))**2 + (Q[:,8]-0.5*(Q[:,2]+Q[:,5]))**2 )**0.5
# rho = (R[:,0]**2 + R[:,1]**2 + R[:,2]**2)**0.5
rho = (R[:,0]**2 + R[:,3]**2 )**0.5

if(too_far):
    print(R[0:3])
    sys.exit("The atoms were too far")  
    

print("min(rho) = ",min(rho))
# plt.figure(3)
# plt.plot(t,rho)

t1 = t/au_s/1E12/1000
dt = []
for i in range(len(t1)-1):
    dt.append(t1[i+1]-t1[i])
print(min(dt))

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

# for i in range(len(t)):
#     E[i][0] = P[i]@P[i]*(au_s*1E12)**2*0.5*29148.94559
# E[:,1] = data[:,-1]

# plt.figure(2)
# plt.plot(t,np.sum(E,axis=-1))
# # plt.plot(t,E[:,0])
# plt.show()

#----------- Radial Acceleration -------------
ar = np.zeros([len(t)-2,4])
# ar = np.zeros([510-10,4])
t_ar = np.zeros((ar.shape[0],))
acc = np.zeros_like(ar)
temp = np.zeros((9,))

rel_acc = np.zeros((3,3))
rel_pos = np.zeros((3,3))
# t = data[:,0]


# t_OP = np.zeros((3,))

# t = data[:,0]

for i in range(len(ar)):
    j = i+1
    hp = t[j+1] -t[j]
    hn = t[j] - t[j-1]
    
    # ar[i][0] = (hn*R[j+1][0] - (hn+hp)*R[j][0] + hp*R[j-1][0]) * 2/(hp*hn*(hp+hn))
    # ar[i][1] = (hn*R[j+1][1] - (hn+hp)*R[j][1] + hp*R[j-1][1]) * 2/(hp*hn*(hp+hn))
    # ar[i][2] = (hn*R[j+1][2] - (hn+hp)*R[j][2] + hp*R[j-1][2]) * 2/(hp*hn*(hp+hn))
    # ar[i][3] = t[j]
    
    for k in range(9):
        temp[k] = (P[j+1][k] - P[j-1][k])/(t[j+1]-t[j-1])
        # temp[k] = (P[j][k] - P[j-1][k])/(t[j]-t[j-1])
    
    # temp = temp/(au_s*1E12)**2
    
    acc[i][0] = (temp[0:3]@temp[0:3])**0.5
    acc[i][1] = (temp[3:6]@temp[3:6])**0.5
    acc[i][2] = (temp[6:9]@temp[6:9])**0.5
    acc[i][3] = t[j]
    
    for k in range(3):
        rel_acc[0] = temp[0:3] - temp[3:6]
        rel_acc[1] = temp[0:3] - temp[6:9]
        rel_acc[2] = temp[3:6] - temp[6:9]
        
        rel_pos[0] = Q[j][0:3] - Q[j][3:6]
        rel_pos[1] = Q[j][0:3] - Q[j][6:9]
        rel_pos[2] = Q[j][3:6] - Q[j][6:9]
        
        # ar[i][k] = -(rel_acc[k]@rel_pos[k])/(rel_pos[k]@rel_pos[k])**0.5
    
    ar[i][3] = t[j]
    
    
    t_ar[i] = t[j]
    
# for ia in range(len(ar)-2):
#     if ( ar[ia+1][0]<-0.5 and ar[ia][0] >=0.5 and ar[ia+2][0] < -0.5):
#         t_OP[0] = t_ar[ia]
#         break
# for ib in range(len(ar)-1):
#     if( ar[ib+1][1]<0.5 and ar[ib][1] >=0.5 and ar[ib+2][1] < -0.5):
#         t_OP[1] = t_ar[ib]
#         break
# for ic in range(len(ar)-1):
#     if( ar[ic+1][2]<0.5 and ar[ic][2] >=0.5 and ar[ic+2][2] < -0.5):
#         t_OP[2] = t_ar[ic]
#         break
# print("t_OP [ps] = ",t_OP)

# print("acc_c = ",acc[ic][2])

t_first = np.zeros((3,))

for j in range(3):
    for i in range(len(acc)):
        if(acc[i][j] > 1.0E-4):
            t_first[j] = acc[i][3]
            break
    if(t_first[j] == 0):
        t_first[j] = acc[-1][3]
print("t_first = ",t_first)

t_2B = min(t_first)
t_3B = max(t_first)
print("t_2B = ",t_2B)
print("t_3B = ",t_3B)

# for i in range(len(acc)):
#     if (acc[i][0] > 5.0):
#         t_2B = acc[i][3]
#         break

# print("t_2B = ",t_2B)

# for i in range(len(acc)):
#     if(acc[i][2] > 5.0):
#         t_3B = acc[i][3]
#         break
    
# # for i in range(10,len(acc)):
# #     dev = 1-acc[i][1]/acc[i][0]
    
# #     if(abs(dev) > 0.1/2 ):
# #         t_3B = acc[i][3]
# #         break
    
# print("t_3B = ",t_3B)
c_id = 0
for i in range(1,len(R)-1):
    if (R[i+1][c_id]>=R[i][0] and R[i-1][c_id]>=R[i][c_id]):
        t_OP = t[i]
        break
print("t_OP = ",t_OP)

ratio = (t_OP-t_3B)/(t_OP-t_2B)
print("ratio = ",ratio)

if(t_3B >= t_OP):
    print(" TSBC ")
elif(ratio < 0.9):
    print("TSBC")
else:
    print("Direct")


# lt = '-'
# plt.figure(2)
# plt.plot(t_ar,ar[:,0],lt,label=r'AB')
# plt.plot(t_ar,ar[:,1],lt,label=r'AC')
# plt.plot(t_ar,ar[:,2],lt,label=r'BC')
# plt.legend()
# plt.ylim([-25,25])    
plt.show()


