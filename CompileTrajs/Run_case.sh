#PBS -q normal 
#PBS -l select=1:ncpus=20:model=ivy
#PBS -N Case_18000
#PBS -o Case_18000.log      
#PBS -l walltime=2:00:00

cd $PBS_O_WORKDIR

Temp=18000
Dir_pref="/nobackupp2/ckondur/CoarseAIR/CoarseAIR_Output/N2O/Recombination/Rates/Pathway/Temp_${Temp}K/T_${Temp}_${Temp}_0_"

case_beg=1
case_end=10

for (( i=$case_beg; i<=$case_end; i++ ))
do
    Dir=$Dir_pref$i/
    echo "Replacing Dir = "$Dir
    sed -i "4s#.*#$Dir#"  Input_files/Files.inp
    
    ./CompileTrajs > Loggers/Log_${Temp}K_$i.log & 
    sleep 5
done

wait 
echo "Done Compiling Trajs for Temp = ${Temp}K"

