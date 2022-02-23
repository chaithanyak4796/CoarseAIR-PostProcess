import numpy as np
import matplotlib.pyplot as plt
import System

sys_name = "N3"
reac_id  = 0

E1 = 100
E2 = 200
Temp = str(E1) + "_" + str(E2)
print(Temp)

#Dir = "/nobackupp2/ckondur/CoarseAIR/CoarseAIR_Output/N3/Recombination/Fixed_omega/" + Temp + "/T_" + Temp + "_"
Dir = "/u/ckondur/CoarseAIR_Output/" + sys_name + "/Recombination/Fixed_omega/" + Temp + "/T_" + Temp + "_"

Omega = np.linspace(-100,100,11,dtype=int)

system = System.System(sys_name)

def Prob_Omega(Dir,Om):
    fname = Dir + str(Om) + "_0/Bins_10_0/trajectories.out"
    Data  = np.loadtxt(fname, skiprows=1)

    b1_tot = Data[:,2]
    b2_tot = Data[:,4]
    v   = Data[:,10]
    j   = Data[:,9]
    arr = Data[:,-3]

    n = len(b1_tot)

    b1 = []
    b2 = []

    for i in range(n):
        #if(v[i] > 0.5 and j[i] > 0.5):
        if((arr[i]-0.5) in system.arr_matrix[reac_id]):
            b1.append(b1_tot[i])
            b2.append(b2_tot[i])

    return len(b1),len(b1)/n

Prob = []
for Om in Omega:
    N,P = Prob_Omega(Dir,Om)
    Prob.append(P)
    #print ("Omega = ",Om/100,"  Prob = ",P)
    print (P)

#plt.figure(1)
#plt.plot(Omega,Prob,label=Temp)
#plt.legend()
#plt.show()
