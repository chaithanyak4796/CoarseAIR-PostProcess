
import numpy as np

Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Fixed_omega/100_200/"
Dir = "/media/chaithanya/Chaithanya_New/QCT_Output/Rate_Constant/Large_Test/Temp_"+str(2000)+"/Statistics/Conv_10/"
# Dir = "/Users/ckondur/Desktop/mount_sshfs/QCT_Output/Rate_Constant/Large_Test/Temp_"+str(2000)+"/Statistics/Conv_40/"

Dir = "/nobackupp2/ckondur/CoarseAIR/CoarseAIR_Output/N3/Recombination/Rates/Temp_10000K/Statistics/"

fname = Dir + "Trajectories-3B-Tot.out"

data = np.loadtxt(fname,skiprows=1,usecols=(-5,-2))#,max_rows=1000)

n_traj = len(data[:,0])

QB_en = []

arr = [17,33,49,19,35,51]

for i in range(n_traj):
    if(data[i][0] in arr):
        QB_en.append(data[i][1])


print("Min QB energy = ",np.min(QB_en))
