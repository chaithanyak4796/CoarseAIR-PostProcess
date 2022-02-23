#PBS -q normal
#PBS -l select=1:ncpus=20:model=ivy
#PBS -N Rates-Path
#PBS -o Rates-Path.log      
#PBS -l walltime=2:00:00

#cd $PBS_O_WORKDIR

#Temp_list=(10000 12000 14000 16000 18000)
Temp_list=(1000 2000 4000 6000 8000)

for Temp in ${Temp_list[@]};
do
    echo "Updating Input files for Temp = "$Temp
    Dir="/u/ckondur/CoarseAIR_Output/N3/Recombination/Rates/Temp_${Temp}K/"
    
    sed -i "4s#.*#$Dir#"   Input_files/Stat_Input.inp
    sed -i "7s#.*#$Temp.0#"  Input_files/Stat_Input.inp
        
    ./Statistics > Loggers/Log-path_${Temp}K.log
    
    echo "Done computing Statistics for Temp = ${Temp}K"
done



