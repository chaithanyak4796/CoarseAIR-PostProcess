import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

Temp = 8000
#Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Test_Current/T_" + str(Temp) + "_" + str(Temp) + "_0_0/"
#fname = Dir + "trajectories-Tot.out"

Dir = "/nobackupp2/ckondur/CoarseAIR/CoarseAIR_Output/N3/Dissociation/Rates/Temp_" + str(Temp) +"K/T_" + str(Temp) + "_" + str(Temp) + "_" + str(1) + "/Bins_1_0/"
fname = Dir + "trajectories.out"

# Read the levels file
levels_fname = './N2_levels.dat'
data_levels = np.loadtxt(levels_fname,skiprows=15,usecols=(0,1,2))
NLevels = len(data_levels)
print(" Number of levels = ",NLevels)

# Read the trajectories-Tot file
data = np.loadtxt(fname,skiprows=1)
conv_int = np.vectorize(int)
v_data = conv_int(np.floor(data[:,5]))
j_data = conv_int(np.floor(data[:,4]))

# Compute the Boltzmann distribution function
kb_Ha = 3.167E-6            # Ha/K
Eint = data_levels[:,2]
beta = 1/(kb_Ha*Temp)

def vj2i(v,j):
    return 500*v+j
def i2vj(i):
    v = i//500
    j = i%500
    return v,j

i_levels = vj2i(data_levels[:,0],data_levels[:,1])
i_inp    = vj2i(v_data,j_data)

Prob_Boltz = (2*data_levels[:,1]+1) * np.exp(-beta*Eint) 
Prob_Boltz = Prob_Boltz/np.sum(Prob_Boltz)
plt.figure(1)
plt.title('Internal Energy Distribution')
plt.semilogy(Eint,Prob_Boltz,'.')

Prob_act = np.zeros((NLevels,))
for n in range(NLevels):
    Prob_act[n] = np.count_nonzero(i_inp == i_levels[n])
Prob_act = Prob_act/sum(Prob_act)
plt.semilogy(Eint,Prob_act,'.')
plt.xlabel('E_int [Ha]')

# Vib_Distrib
n_vib = len(np.unique(data_levels[:,0]))
Prob_Boltz_vib = np.zeros((n_vib,))
E_vib = np.zeros_like(Prob_Boltz_vib)
for i in range(NLevels):
    v = int(data_levels[i][0])
    j = int(data_levels[i][1])
    Prob_Boltz_vib[v] += Prob_Boltz[i]
    if(j == 0):
        E_vib[v] = Eint[i]
    
plt.figure(2)
plt.title('Vib Energy Distribution')
plt.semilogy(E_vib,Prob_Boltz_vib)

Prob_act_vib = np.zeros((n_vib,))
for v in range(n_vib):
    Prob_act_vib[v] = np.count_nonzero(v_data==v)
Prob_act_vib = Prob_act_vib/np.sum(Prob_act_vib)
plt.semilogy(E_vib,Prob_act_vib)
plt.xlabel('E_vib [Ha]')
    
plt.show()
