import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy.linalg as la
import seaborn as sns

E1 = 100
E2 = 200
Om = 30

Temp = str(E1) + "_" + str(E2)

fname = "/nobackupp2/ckondur/CoarseAIR/CoarseAIR_Output/N3/Recombination/Fixed_omega/" + Temp + "/T_" + Temp + "_" + str(Om) + "_0/Bins_10_0/trajectories.out"

if(1):
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['text.usetex'] = 'true'
    plt.rcParams['font.size'] = 20
    plt.rcParams['axes.labelsize'] = 15
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titlesize'] = 30
    plt.rcParams['xtick.labelsize'] = 15
    plt.rcParams['ytick.labelsize'] = 15
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['lines.markersize'] = 5
    plt.rcParams['axes.grid'] = False
    plt.rcParams['axes.grid.which'] = 'both'
    plt.rcParams['grid.linestyle'] = '--'
    plt.rcParams['grid.linewidth'] = 1

    plt.rc('font',weight='bold')
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

Data = np.loadtxt(fname, skiprows=1)

b1_max = np.unique(Data[:,2])
b2_max = np.unique(Data[:,4])

#b1_max = np.insert(b1_max,0,0)
#b2_max = np.insert(b2_max,0,0)

idx_rec = np.where(Data[:,9]!=0.5)[0]

b1 = Data[idx_rec,3]
b2 = Data[idx_rec,5]

n_recomb = len(b1)
print("n_recomb = ",n_recomb)

n_traj = len(Data)
print("n_traj   = ",n_traj)

hist_rec,xedg,yedg = np.histogram2d(b1, b2, bins=[b1_max,b2_max])
hist_tot,xedg,yedg = np.histogram2d(Data[:,2]*0.95, Data[:,4]*0.95, bins=[b1_max,b2_max])

Prob = hist_rec/hist_tot
#print(Prob)

temp1 = b1_max[:-1]
temp2 = b2_max[:-1]
xx,yy = np.meshgrid(temp1,temp2)

x,y = xx.ravel(),yy.ravel()
top = Prob.ravel()

#print(len(x),len(y),len(top))

bottom = np.zeros_like(top)
width = depth = 1
fig = plt.figure(figsize=(8,6))
ax = Axes3D(fig)

ax.bar3d(y, x, bottom, width, depth, top, shade=True,color=None)
ax.view_init(azim=45,elev=20)
ax.set_xlabel(r"$b_1 [a_0] $ ")
ax.set_ylabel(r"$b_2 [a_0] $")
#ax.set_xticks(b1_max)
#ax.set_yticks(b2_max)

plt.savefig('./Prob_imp_'+str(Om)+'.eps',forat='eps')

plt.show()
