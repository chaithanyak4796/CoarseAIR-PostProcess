#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 10:01:31 2019

@author: chaithanya
"""

import numpy as np
import matplotlib.pyplot as plt
# from mpldatacursor import datacursor

if(0):
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['text.usetex'] = 'true'
    plt.rcParams['font.size'] = 15
    plt.rcParams['axes.labelsize'] = 10
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titlesize'] = 35
    plt.rcParams['xtick.labelsize'] = 15
    plt.rcParams['ytick.labelsize'] = 15
    plt.rcParams['lines.linewidth'] = 1
    plt.rcParams['lines.markersize'] = 10
    plt.rcParams['axes.grid'] = False
    plt.rcParams['axes.grid.which'] = 'both'
    plt.rcParams['grid.linestyle'] = '--'
    plt.rcParams['grid.linewidth'] = 1

#plt.rc('font',weight='bold')
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

# plt.close('all')
case = 2  # 0 : QCT , 1 : Stat

if(case == 1):
    
    # Dir = "/Users/ckondur/Desktop/mount_sshfs/QCT_Output/Rate_Constant/Large_Test/Temp_" + str(6000) +"/Statistics/Conv_10/"
    Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Rate_Constant/Large_Test/Temp_" + str(6000) +"/Statistics/Conv_10/"
    # Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Fixed_omega/100_100/T_100_100_60_0/Statistics/"
    fname  = Dir + "Trajectories_tot.out"
    
    Data = np.loadtxt(fname, skiprows=1)

    b1_tot = Data[:,6]
    b2_tot = Data[:,8]
    v  = Data[:,11]
    j  = Data[:,12] 
    
    n = len(b1_tot)
    # b1 = b1_tot
    # b2 = b2_tot

    b1 = []
    b2 = []
    
    for i in range(n):
        if(v[i] >=0 and j[i] >= 0):
            b1.append(b1_tot[i])
            b2.append(b2_tot[i])

elif(case==0):
    #Dir = "/scratch1/07819/ckondur/CG-QCT/Omega_Fixed/Prob/100_100/T_100_100_0_0/"
    Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Test_Current/T_100_200_60_0/"
    #Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Testing/N2/T_100_200_0_0/"
    fname = Dir + "trajectories.out"
    
    Data = np.loadtxt(fname, skiprows=1)

    b1_tot = Data[:,2]
    b2_tot = Data[:,4]
    v  = Data[:,10]
    j  = Data[:,9] 

    n = len(b1_tot)
    
    b1 = []
    b2 = []

    for i in range(n):
        if(v[i] > 0.5 and j[i] > 0.5):
            b1.append(b1_tot[i])
            b2.append(b2_tot[i])
    
    print("n_recomb = ",len(b1))
    print("P_recomb = ",len(b1)/n)

def Prob_Omega(Dir,Om):
    fname = Dir + str(Om) + "_0/Bins_10_0/trajectories.out"
    Data = np.loadtxt(fname, skiprows=1)

    b1_tot = Data[:,2]
    b2_tot = Data[:,4]
    v      = Data[:,10]
    j      = Data[:,9]
    arr    = Data[:,11]

    n = len(b1_tot)

    b1 = []
    b2 = []

    reac_arr = np.array([[16.5,17.5,19.5],[32.5,33.5,35.5],[48.5,49.5,51.5]])

    for i in range(n):
        #if(v[i] > 0.5 and j[i] > 0.5):
        if(arr[i] in reac_arr[0]):
            b1.append(b1_tot[i])
            b2.append(b2_tot[i])
            
    return len(b1),len(b1)/n

E1 = 100
E2 = 200
Temp = str(E1) + "_" + str(E2)
#Dir="/media/chaithanya/Chaithanya_New/QCT_Output/Test_Current/T_"+Temp+"_"
#Dir = "/nobackupp2/ckondur/CoarseAIR/CoarseAIR_Output/N2O/Recombination/Fixed_omega/" + Temp + "/T_" + Temp + "_"
#Dir = "/nobackupp2/ckondur/CoarseAIR/CoarseAIR_Output/N2O/Recombination/Fixed_omega/" + Temp + "/T_" + Temp + "_"

Dir = "/u/ckondur/CoarseAIR_Output/N2O/Recombination/Fixed_omega/" + Temp + "/T_" + Temp + "_"

#Omega = [-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100]
Omega = [-100, -80, -60, -40, -20, -10, 0, 10, 20, 40, 60, 80, 100]

Prob = []
if(case == 2):
    for Om in Omega:
        N,P = Prob_Omega(Dir,Om)
        Prob.append(P)

        print ("Omega = ",Om/100,"  Prob = ",P)

    plt.figure()
    plt.plot(Omega,Prob,'-o')
    plt.show()

plot = 0
if(plot):
    bin_arr = np.linspace(0,20,9)
    fig = plt.figure(figsize=(10,8))
    plt.hist2d(b1,b2,density=True, bins=[bin_arr,bin_arr])
    plt.ylim((0,16))
    plt.xlim((0,16))
    # plt.title(r'$d_1 = 80 Bo $')
    plt.xlabel(r'$b_1 [Bo] $')
    plt.ylabel(r'$b_2 [Bo] $')
    #plt.hist2d(b1_tot,b2_tot,bins=8)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=25) 
    # datacursor(draggable=True)

    locs, labels = plt.xticks()
    plt.xticks(bin_arr)

    locs, labels = plt.yticks()
    plt.yticks(bin_arr)

    plt.show()
