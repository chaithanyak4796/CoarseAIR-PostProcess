Temp_list=(10000 12000 14000 16000 18000)

for T in ${Temp_list[@]};
do
    tail -n+2 Temp_${T}K/Statistics/Total_Rate_Constant_1.out >> Poission_yes_1.out
    tail -n+2 Temp_${T}K/Statistics/Total_Rate_Constant_2.out >> Poission_yes_2.out
done
