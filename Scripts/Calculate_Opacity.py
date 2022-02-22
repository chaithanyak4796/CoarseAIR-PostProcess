#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 15:47:32 2022

@author: chaithanya
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
import System

plt.close('all')

Temp         = 18000
sys_name     = "N2O"
reac_id      = 0
resolve_path = True
num_bins     = 20

Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Results/N2O/Recombination/Rates/Pathway/Temp_" + str(Temp) + "K/Statistics/"
# Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Results/N3/Recombination/Rates/Temp_" + str(Temp) + "K/Statistics/"

fname = Dir + "Trajectories-3B-Tot.out"

Data   = np.loadtxt(fname,skiprows=1)
system = System.System(sys_name)

b1 = Data[:,6]
b2 = Data[:,8]

b1_max = Data[:,5]
b2_max = Data[:,7]

arr = Data[:,11]

rings_1 = np.linspace(0, max(b1_max), num_bins)
rings_2 = np.linspace(0, max(b2_max), num_bins)

hist_tot_b1, edg_b1 = np.histogram(b1, bins=rings_1)
hist_tot_b2, edg_b2 = np.histogram(b2, bins=rings_2)

if(0):
    idx_reac = np.where(Data[:,10] > 0)[0]
    if(system.num_reac > 1):
        idx_reac_t = []
        for i in range(len(system.arr_matrix[reac_id])):
            temp = np.where(arr[idx_reac] == system.arr_matrix[reac_id][i])[0]
            idx_reac_t.append(idx_reac[temp])

        idx_reac = np.concatenate(idx_reac_t,axis=0)
else:
     idx_reac = []
     if(resolve_path):
         idx_reac_path = [ [], [], [] ]
     for i in range(len(Data)):
         if (arr[i] in system.arr_matrix[reac_id]):
             idx_reac.append(i)

             if(resolve_path):
                 idx_reac_path[int(Data[i][13])-1].append(i)

     idx_reac = np.array(idx_reac)
     if(resolve_path):
         for i in range(len(idx_reac_path)):
             idx_reac_path[i] = np.array(idx_reac_path[i])

hist_reac_b1 , edg_b1 = np.histogram(b1[idx_reac], bins = edg_b1)
hist_reac_b2 , edg_b2 = np.histogram(b2[idx_reac], bins = edg_b2)

Op_b1 = hist_reac_b1/hist_tot_b1
Op_b2 = hist_reac_b2/hist_tot_b2

if(resolve_path):
    hist_reac_b1_path_L, edg_b1 = np.histogram(b1[idx_reac_path[0]], bins = edg_b1)
    hist_reac_b2_path_L, edg_b2 = np.histogram(b2[idx_reac_path[0]], bins = edg_b2)

    hist_reac_b1_path_C, edg_b1 = np.histogram(b1[idx_reac_path[1]], bins = edg_b1)
    hist_reac_b2_path_C, edg_b2 = np.histogram(b2[idx_reac_path[1]], bins = edg_b2)

    hist_reac_b1_path_D, edg_b1 = np.histogram(b1[idx_reac_path[2]], bins = edg_b1)
    hist_reac_b2_path_D, edg_b2 = np.histogram(b2[idx_reac_path[2]], bins = edg_b2)

    Op_b1_L = hist_reac_b1_path_L/hist_tot_b1
    Op_b2_L = hist_reac_b2_path_L/hist_tot_b2
    Op_b1_C = hist_reac_b1_path_C/hist_tot_b1
    Op_b2_C = hist_reac_b2_path_C/hist_tot_b2
    Op_b1_D = hist_reac_b1_path_D/hist_tot_b1
    Op_b2_D = hist_reac_b2_path_D/hist_tot_b2

    # Scaling the Opacity functions
    Op_b1_L /= np.sum(Op_b1)
    Op_b2_L /= np.sum(Op_b2)
    Op_b1_C /= np.sum(Op_b1)
    Op_b2_C /= np.sum(Op_b2)
    Op_b1_D /= np.sum(Op_b1)
    Op_b2_D /= np.sum(Op_b2)
    
rings_b1 = np.zeros(len(Op_b1))
rings_b2 = np.zeros(len(Op_b2))

# Scaling the Opacity functions
Op_b1 /= np.sum(Op_b1)
Op_b2 /= np.sum(Op_b2)

for i in range(len(rings_b1)):
    rings_b1[i] = (edg_b1[i] + edg_b1[i+1])/2
    
for i in range(len(rings_b2)):
    rings_b2[i] = (edg_b2[i] + edg_b2[i+1])/2
    


if not resolve_path:
    plt.figure(figsize=(12,8))
    plt.subplot(121)
    plt.plot(rings_b1, Op_b1)
    plt.xlabel('b1 [A]')
    plt.ylabel(r'$P_r(b1)/P_r$')
    
    plt.subplot(122)
    plt.plot(rings_b2, Op_b2)
    plt.xlabel('b2 [A]')
    plt.ylabel(r'$P_r(b2)/P_r$')
    prefix  = Dir + "Opacity_reac_" + str(reac_id+1) + "_"

    fname_1 = prefix + str(1) + ".dat"
    fname_2 = prefix + str(2) + ".dat"
    
    fw1 = open(fname_1, "w")
    fw2 = open(fname_2, "w")
    
    #fw1.write("   b1 [A]       Op_b1\n")
    #fw2.write("   b2 [A]       Op_b2\n")
    
    fw1.write("Variables = \"b<sub>1<\sub>\", \"P<sub>r</sub>(b<sub>1</sub>)/P<sub>r</sub>\" \n")
    fw1.write("Zone T = \"%dK\", I = %d, J = 1, K= 1, F = POINT, DT = (DOUBLE, DOUBLE) \n"%(Temp, len(rings_b1)))
    
    fw2.write("Variables = \"b<sub>2<\sub>\", \"P<sub>r</sub>(b<sub>2</sub>)/P<sub>r</sub>\" \n")
    fw2.write("Zone T = \"%dK\", I = %d, J = 1, K= 1, F = POINT, DT = (DOUBLE, DOUBLE) \n"%(Temp, len(rings_b2)))
    
    for i in range(len(rings_b1)):
        fw1.write("%8.6E  %8.6E\n"%(rings_b1[i], Op_b1[i]))
    for i in range(len(rings_b2)):
        fw2.write("%8.6E  %8.6E\n"%(rings_b2[i], Op_b2[i]))
    
    fw1.close()
    fw2.close()

if(resolve_path):
    plt.figure(2,figsize=(12,8))
    plt.subplot(121)
    plt.plot(rings_b1, Op_b1_L, label="L")
    plt.plot(rings_b1, Op_b1_C, label="C")
    plt.plot(rings_b1, Op_b1_D, label="D")
    plt.xlabel('b1 [A]')
    plt.ylabel(r'$P_r(b1)/P_r$')
    plt.legend()

    plt.subplot(122)
    plt.plot(rings_b2, Op_b2_L, label="L")
    plt.plot(rings_b2, Op_b2_C, label="C")
    plt.plot(rings_b2, Op_b2_D, label="D")
    plt.xlabel('b2 [A]')
    plt.ylabel(r'$P_r(b2)/P_r$')
    plt.legend()

    prefix  = Dir + "Opacity_path_reac_" + str(reac_id+1) + "_"

    fname_1 = prefix + str(1) + ".dat"
    fname_2 = prefix + str(2) + ".dat"

    fw1 = open(fname_1, "w")
    fw2 = open(fname_2, "w")

    #fw1.write("   b1 [A]       Op_b1_L	      Op_b1_C  	    Op_b1_D    	  Op_b1_T\n")
    #fw2.write("   b2 [A]       Op_b2_L	      Op_b2_C  	    Op_b2_D    	  Op_b2_T\n")
    
    fw1.write("Variables = \"b<sub>1<\sub>\", \"Lindemann\", \"Chaperon\", \"Direct\", \"Total\" \n")
    fw1.write("Zone T = \"%dK\", I = %d, J = 1, K= 1, F = POINT, DT = (DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE) \n"%(Temp, len(rings_b1)))
    
    fw2.write("Variables = \"b<sub>2<\sub>\", \"Lindemann\", \"Chaperon\", \"Direct\", \"Total\" \n")
    fw2.write("Zone T = \"%dK\", I = %d, J = 1, K= 1, F = POINT, DT = (DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE) \n"%(Temp, len(rings_b2)))
    
    for i in range(len(rings_b1)):
        fw1.write("%8.6E  %8.6E  %8.6E  %8.6E  %8.6E\n"%(rings_b1[i], Op_b1_L[i], Op_b1_C[i], Op_b1_D[i], Op_b1_L[i]+Op_b1_C[i]+Op_b1_D[i]))
    for i in range(len(rings_b2)):
        fw2.write("%8.6E  %8.6E  %8.6E  %8.6E  %8.6E\n"%(rings_b2[i], Op_b2_L[i], Op_b2_C[i], Op_b2_D[i], Op_b2_L[i]+Op_b2_C[i]+Op_b2_D[i]))
    fw1.close()
    fw2.close()
    


plt.show()
