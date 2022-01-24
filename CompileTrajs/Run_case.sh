Temp=12000
Dir_pref="/media/chaithanya/Chaithanya_New/QCT_Output/Results/N3/Recombination/Rates/Pathway/Temp_${Temp}K/T_${Temp}_${Temp}_0_"

case_beg=1
case_end=10

for (( i=$case_beg; i<=$case_end; i++ ))
do
    Dir=$Dir_pref$i/
    echo "Replacing Dir = "$Dir
    sed -i "4s#.*#$Dir#"  Input_files/Files.inp
    
    ./CompileTrajs > Loggers/Log-new_${Temp}K_$i.log & 
    sleep 5
done

wait 
echo "Done Compiling Trajs for Temp = ${Temp}K"

