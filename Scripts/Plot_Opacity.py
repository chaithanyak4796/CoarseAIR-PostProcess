import numpy as np
import matplotlib.pyplot as plt

#Temp_list = [10000]
#Temp_list = [10000,12000,14000,16000,18000]
Temp_list=[1000,2000,4000,8000,10000,14000,18000]

for i in range(len(Temp_list)):
    Temp = Temp_list[i]
    Dir="/media/chaithanya/Chaithanya_New/QCT_Output/Results/N3/Recombination/Rates/Temp_" + str(Temp) + "K/Statistics/"

    fname1 = Dir + "Opacity_Func1.out"
    fname2 = Dir + "Opacity_Func2.out"

    OP1 = np.loadtxt(fname1,skiprows=1)
    OP2 = np.loadtxt(fname2,skiprows=1)

    b1 = OP1[:,0]
    b2 = OP2[:,0]

    OP1 = OP1[:,3]/np.sum(OP1[:,3])
    OP2 = OP2[:,3]/np.sum(OP2[:,3])

    if(i==0):
        Data_1 = [OP1]
        Data_2 = [OP2]
    else:
        Data_1.append(OP1)
        Data_2.append(OP2)
        
    plt.figure(1,figsize=(16,12))
    plt.subplot(121)
    plt.plot(b1, OP1, '-o', label = str(Temp)+"K")
    plt.xlabel(r"$b_1 [Bo]$")
    plt.ylabel(r"$P(b_1|recomb)/P(recomb)$")
    plt.legend()
    
    plt.subplot(122)
    plt.plot(b2, OP2, '-o', label = str(Temp)+"K")
    plt.xlabel(r"$b_2 [Bo]$")
    plt.ylabel(r"$P(b_2|recomb)/P(recomb)$")
    plt.legend()


#print(Data_2)
#print(OP2)

fname1 = 'Opacity_1.dat'
fname2 = 'Opacity_2.dat'

f1 = open(fname1,"w")
for i in range(len(b1)):
    f1.write("%6.3f  "%(b1[i]))
    for j in range(len(Temp_list)):
        f1.write("%6.4E  "%(Data_1[j][i]))
    f1.write("\n")

f2 = open(fname2,"w")
for i in range(len(b2)):
    f2.write("%6.3f  "%(b2[i]))
    for j in range(len(Temp_list)):
        f2.write("%6.4E  "%(Data_2[j][i]))
    f2.write("\n")

f1.close()
f2.close()


plt.show()
