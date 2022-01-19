import numpy as np
import sys

system = "N3"

Kb = 3.166829E-6  # Ha/K
if (system == "O3"):
    m1 = 29148.94559  # me
elif(system == "N3"):
    m1 = 25526.04298

m2 = m1
m3 = m1

au_s = 2.4188843265E-17
bo_Angs = 0.529177249
Angs_m = 1E-10

if (len(sys.argv) == 2):
    Temp = int(sys.argv[1])
else:
    Temp = 10000

case = [1,2,3,4,5,6,7,8,9,10]

for c in range(len(case)):
    # Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Test_Current/T_" + str(Temp) + "_" + str(Temp) + "_0_" + str(case[c]) + "/"
    # Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Testing/Dissociation/T_" + str(Temp) + "_" + str(Temp) + "_0_0/"
    #Dir = "/nobackupp2/ckondur/CG_QCT/CG_QCT_Output/N2/Dissociation/Temp_" + str(Temp) +"K/T_" + str(Temp) + "_" + str(Temp) + "_0_" + str(case[c]) + "/"
    #fname = Dir + "trajectories-Tot.out"

    Dir = "/nobackupp2/ckondur/CoarseAIR/CoarseAIR_Output/" + system+ "/Dissociation/Rates/Temp_" + str(Temp) +"K/T_" + str(Temp) + "_" + str(Temp) + "_" + str(case[c]) + "/Bins_1_0/"
    fname = Dir + "trajectories.out"

    data_c = np.loadtxt(fname,skiprows=1,usecols=(0,2,3,-3))
    
    if(c==0):
        data = data_c
    else:
        data = np.concatenate((data,data_c))

arr = data[:,3]

diss = [18.50, 34.50, 50.50]

idx1 = np.where(arr == diss[0])[0]
idx2 = np.where(arr == diss[1])[0]
idx3 = np.where(arr == diss[2])[0]

num_diss = len(idx1) + len(idx2) + len(idx3)

print("Number of diss = ",num_diss)

b_max = data[:,1]
b     = data[:,2]

rings   = np.unique(b_max)
n_rings = len(rings)
rings   = np.insert(rings,0,0)
#print(" Number of rings = ",n_rings)
#print(" Rings : ",rings)

Area = np.zeros(n_rings,)
for i in range(len(Area)):
    Area[i] = np.pi*(rings[i+1]**2 - rings[i]**2)
    
diss_idx = np.concatenate((idx1,idx2,idx3))
b_diss   = data[diss_idx][:,2]

hist_tot, edg = np.histogram(b,bins=rings)
hist_diss,edg = np.histogram(b_diss,bins=rings)
#print(hist_tot)#/np.sum(hist_tot))
# print(edg)

# count_all = np.zeros(n_rings)
# for i in range(data.shape[0]):
#     for j in range(n_rings):
#         if(b[i]>=rings[j] and b[i]<=rings[j+1]):
#             count_all[j] += 1

# print(count_all)

Prob = hist_diss/hist_tot

k     = 0
k_std = 0

mu = (m1+m2)*m3/(m1+m2+m3)

for i in range(n_rings-1):
    k += Area[i]*Prob[i]
    
    temp   = Prob[i]*(1-Prob[i])/hist_tot[i]
    k_std += Area[i]*temp**0.5


coeff =  (8*Kb*Temp/np.pi/mu)**0.5
conv_fac = (bo_Angs*Angs_m*100)**3/au_s

if(system == "O3"): 
    adiabatic_corr = 16/3
    symm_fac = 1.0

else: 
    adiabatic_corr = 1.0
    symm_fac = 1.0

k     *= coeff * conv_fac * adiabatic_corr * symm_fac
k_std *= coeff * conv_fac * adiabatic_corr * symm_fac

print("k_diss [cm^3/s] = %6.8E"%(k))
print("k_std  [cm^3/s] = %6.8E"%(k_std))
