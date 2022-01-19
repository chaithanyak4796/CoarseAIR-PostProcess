import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from sklearn.neighbors import KernelDensity

case = 1 # 0: Total; 1: Path_specific

Temp_list = [4000, 16000]
E_shift = -19.83475060613458 + 9.917

def kde_fit(bw,x,y):
    kde = KernelDensity(bandwidth=bw,kernel='gaussian')
    kde.fit(y[:,None])
    logprob = kde.score_samples(x[:,None])
    return np.exp(logprob)
                
for i in range(len(Temp_list)):
    Temp = Temp_list[i]
    Dir="/media/chaithanya/Chaithanya_New/QCT_Output/Results/N3/Recombination/Rates/Temp_" + str(Temp) + "K/Statistics/"
    #Dir="/media/chaithanya/Chaithanya_New/QCT_Output/Test_Current/Statistics/"
    #Dir="/media/chaithanya/Chaithanya_New/QCT_Output/testing/CoarseAIR/Rates/Temp_" + str(Temp) + "K/Statistics/"
    

    if(case == 0):
        fname = Dir + "Trajectories-3B-Tot.out"
        data    = np.loadtxt(fname,skiprows=1,usecols=(-3,-2))
        idx_rec = np.where(data[:,0]>0)[0]

        E_int_data = data[idx_rec,1] + E_shift
        nrec       = len(E_int_data)

        print("\nTemp          = ",Temp)
        print("Number of rec = ",nrec)
        print("Minimum E_int = ", min(E_int_data))
        save_name = Dir + "Internal_Energy_Dist.dat"

    if(case == 1):
        Dir="/media/chaithanya/Chaithanya_New/QCT_Output/Results/N3/Recombination/Rates/Pathway/Temp_" + str(Temp) + "K/Statistics_1mil/"
        fname = Dir + "Trajectories_tot.out"
        data    = np.loadtxt(fname,skiprows=1,usecols=(-4,-3))
        L_data  = data[np.where(data[:,0]==1)[0],1] + E_shift
        C_data  = data[np.where(data[:,0]==2)[0],1] + E_shift
        D_data  = data[np.where(data[:,0]==3)[0],1] + E_shift
        E_int_data  = np.concatenate((L_data,C_data,D_data))
        nrec    = len(L_data) + len(C_data) + len(D_data)

        save_name = Dir + "Internal_Energy_Dist_path.dat"

        print("Number of L = ",len(L_data))
        print("Number of C = ",len(C_data))
        print("Number of D = ",len(D_data))
        print("Total number of rec = ",len(E_int_data))
        
    nbins = 100
    en = np.linspace(0,max(E_int_data),nbins)
    bin_wid = en[1] - en[0]

    if(case == 0):
        kde = KernelDensity(bandwidth=bin_wid,kernel='gaussian')
        kde.fit(E_int_data[:,None])
        logprob = kde.score_samples(en[:,None])
        kde_data = np.array([en,np.exp(logprob)]).T

        ax_tot = plt.semilogy(en,np.exp(logprob))
        #ax_tot.set_yscale('log')
        
        head = "E_int[eV]   f(E_int)"
        np.savetxt(save_name,kde_data,fmt='%6.4E',delimiter='  ',header=head)

    else:
        kde_L = kde_fit(bin_wid,en,L_data)
        kde_C = kde_fit(bin_wid,en,C_data)
        kde_D = kde_fit(bin_wid,en,D_data)
        kde_T = kde_fit(bin_wid,en,E_int_data)

        plt.semilogy(en,kde_L)
        plt.semilogy(en,kde_C)
        plt.semilogy(en,kde_D)

        kde_data = np.array([en,kde_L,kde_C,kde_D,kde_T]).T
        head = "E_int[eV]   f_L(E_int)  f_C(E_int)  f_D(E_int)  f_T(E_int)"
        np.savetxt(save_name,kde_data,fmt='%6.4E',delimiter='  ',header=head)
        
plt.ylim([1E-4,10.0])
plt.show()
