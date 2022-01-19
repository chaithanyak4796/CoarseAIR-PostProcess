#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 11:20:29 2021

@author: chaithanya
"""
import numpy as np

Temp   = 8000
NProcs = 1 

# Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Test_Current/T_" + str(Temp) + "_" + str(Temp) + "_0_0/Bins_1951_0/Node_1/Proc_"
Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Testing/Dissociation/T_" + str(Temp) + "_" + str(Temp) + "_0_0/Bins_1_0/Node_1/Proc_"

Dir = "/nobackupp2/ckondur/CG_QCT/CG_QCT_Output/Test_Current/T_" + str(Temp) + "_" + str(Temp) + "_0_0/Bins_1_0/Node_1/Proc_"

iP = 1
iT = 204

fname = Dir + str(iP) + "/PaQSol.out"
data  = np.loadtxt(fname,skiprows=1)

idx  = np.where(data[:,0] == iT)[0][0]
PaQ  = data[idx][3:15]

P = np.zeros(9)
Q = np.zeros(9)

P[0:6] = PaQ[0:6]
Q[0:6] = PaQ[6:12]
for i in range(3):
    P[i+6] = -(P[i] + P[i+3])
    Q[i+6] = -(Q[i] + Q[i+3])

QA = (Q[0:3]+Q[3:6])/2
QB = Q[6:9]
PA = (P[0:3]+P[3:6])/2
PB = P[6:9]

R = QB-QA
V = PB-PA

ct = (R@V)**2/( ( R@R) * (V@V) )
ct = min(1.0,ct)

b1 = (R@R)**0.5 * (1-ct)**0.5
print(b1)
