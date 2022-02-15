import numpy as np

Temp = np.array([1000, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000])

pref     = "Poisson_yes_"
num_reac = 2
pathways = False

T_rate = np.zeros((num_reac,len(Temp),2))
L_rate = np.zeros_like(T_rate)
C_rate = np.zeros_like(T_rate)
D_rate = np.zeros_like(T_rate)

for i in range(len(Temp)):
    T = Temp[i]
    
    for j in range(num_reac):
        fname = "./Temp_" + str(T) + "K/Statistics/Total_Rate_Constant_" + str(j+1) + ".out"
        
        data = np.loadtxt(fname,skiprows=1,usecols=(2,3))
    
        if(pathways):
            L_rate[j][i][0] = data[0][0]
            L_rate[j][i][1] = data[0][1]
            C_rate[j][i][0] = data[1][0]
            C_rate[j][i][1] = data[1][1]
            D_rate[j][i][0] = data[2][0]
            D_rate[j][i][1] = data[2][1]
            T_rate[j][i][0] = data[3][0]
            T_rate[j][i][1] = data[3][1]
        else:
            T_rate[j][i][0] = data[0]
            T_rate[j][i][1] = data[1]

        #print(T_rate[j])

if(pathways):
    header = "Temp       Rate_L         Std_L         Rate_C         Std_C         Rate_D        Std_D         Rate_Tot       Std_tot\n"
    for i in range(num_reac):
        fname = pref + str(i+1) + ".out"
        fw    = open(fname,"w")
        fw.write(header)
        for j in range(len(Temp)):
            fw.write("%.2f  %12.6E  %12.6E  %12.6E  %12.6E  %12.6E  %12.6E  %12.6E  %12.6E\n"%(Temp[j], L_rate[i][j][0], L_rate[i][j][1], C_rate[i][j][0], C_rate[i][j][1],
                                                                                       D_rate[i][j][0], D_rate[i][j][1],  T_rate[i][j][0], T_rate[i][j][1]))
    fw.close()

else:
    header = "Temp       Rate_Tot       Std_tot\n"
    for i in range(num_reac):
        fname = pref + str(i+1) + ".out"
        fw    = open(fname,"w")
        fw.write(header)
        for j in range(len(Temp)):
            fw.write("%.2f  %12.6E  %12.6E\n"%(Temp[j], T_rate[i][j][0], T_rate[i][j][1])) 
