import numpy as np
import matplotlib.pyplot as plt

Temp_list = [1000,2000,4000,6000,8000,10000,12000,14000,16000,18000]

for i in range(len(Temp_list)):
    Temp = Temp_list[i]
    Dir="/media/chaithanya/Chaithanya_New/QCT_Output/Results/N3/Recombination/Rates/Temp_" + str(Temp) + "K/Statistics/"
    #Dir="/media/chaithanya/Chaithanya_New/QCT_Output/Test_Current/Statistics/"
    #Dir="/media/chaithanya/Chaithanya_New/QCT_Output/testing/CoarseAIR/Rates/Temp_" + str(Temp) + "K/Statistics/"
    fname = Dir + "FinStates.out"

    Data = np.loadtxt(fname,skiprows=1)

    E   = Data[:,5]
    n_i = Data[:,3]
    k_i = Data[:,7]
    v   = Data[:,1]
    j   = Data[:,2]

    N = int(np.sum(n_i))
    print("Total number of states = ",len(k_i))

    v_unq = np.unique(v)
    j_unq = np.unique(j)

    k_v = np.zeros_like(v_unq)
    k_j = np.zeros_like(j_unq)
    N_v = np.zeros_like(v_unq)
    N_j = np.zeros_like(j_unq)


    for i in range(len(k_v)):
        idx = np.where(v == v_unq[i])[0]
        k_v[i] = np.sum(k_i[idx])
        N_v[i] = np.sum(n_i[idx])/N

    for i in range(len(k_j)):
        idx = np.where(j == j_unq[i])[0]
        k_j[i] = np.sum(k_i[idx])
        N_j[i] = np.sum(n_i[idx])/N

    v_data = np.array([v_unq,N_v]).T
    j_data = np.array([j_unq,N_j]).T

    fname_v = Dir + "/v_distrib.dat"
    fname_j = Dir + "/j_distrib.dat"

    np.savetxt(fname_v,v_data,fmt='%.8f')
    np.savetxt(fname_j,j_data,fmt='%.8f')

    plt.subplot(121)
    plt.semilogy(v_unq,N_v,label=str(Temp))

    plt.subplot(122)
    plt.semilogy(j_unq,N_j,label=str(Temp))

plt.legend()
plt.show()

