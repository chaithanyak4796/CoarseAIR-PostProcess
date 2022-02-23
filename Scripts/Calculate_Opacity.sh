#PBS -q ldan
#PBS -l select=1:ncpus=20
#PBS -N Opacity
#PBS -o Opacity.log      
#PBS -l walltime=4:00:00

#cd $PBS_O_WORKDIR

#Temp_list=(10000 12000 14000 16000 18000 1000 2000 4000 6000 8000)
Temp_list=(10000 12000 14000 16000 18000)

module load comp-intel/2018.3.222 python3/Intel_Python_3.6_2018.3.222

for Temp in ${Temp_list[@]};
do
    echo "Temp = "$Temp
    python Calculate_Opacity.py $Temp
done

#wait
