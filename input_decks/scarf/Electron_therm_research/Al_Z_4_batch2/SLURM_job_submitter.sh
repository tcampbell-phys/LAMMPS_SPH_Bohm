#!/bin/bash

### SLURM job submissions

start_string="in.placeholder"
start_deck_var="variable out_number string "
echo "start deck letter:"
read start_deck_let
base_input="in.scarf_e_therm_Z_4_gamma_scan_Nepe_16_"

echo "Number range start: "
read start

echo "Number range end: "
read end

echo "Letter variable (0 = a, 1 = b, 2 = c, 3 = d): "
read letter_var

if [ $letter_var -eq 0 ];
then
    suffix="a"
fi

if [ $letter_var -eq 1 ];
then
    suffix="b"
fi

if [ $letter_var -eq 2 ];
then
    suffix="c"
fi

if [ $letter_var -eq 3 ];
then
    suffix="d"
fi

echo "selected suffix: ${suffix}"


for ((i=$start; i<=$end; i++))
do
    if [ $i -eq ${start} ];
    then
    echo "Start run_lammps_SPH_Bohm string: ${start_string}"
    echo "Replacement run_lammps_SPH_Bohm string: ${base_input}${i}${suffix}"
    sed -i -e "s/< ${start_string}/< ${base_input}${i}${suffix}/g" ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/ele_therm_research/run_lammps_SPH_Bohm
    
    echo "Start in.* string: ${start_deck_var}${i}${start_deck_let}"
    echo "Replacement in.* string: ${start_deck_var}${i}${suffix}"
    cp ${base_input}${i}a ${base_input}${i}${suffix}
    sed -i -e "s/${start_deck_var}${i}${start_deck_let}/${start_deck_var}${i}${suffix}/g" ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/ele_therm_research/${base_input}${i}${suffix}

    sbatch run_lammps_SPH_Bohm
    final_out=${base_input}${i}${start_deck_let}
    #cp ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/Al_research/3_5eV_5_2gcc/run_lammps_SPH_Bohm ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/Al_research/3_5eV_5_2gcc/run_lammps_SPH_Bohm${i}
    fi
    if [ $i -ne ${start} ];
    then
    echo "Start run_lammps_SPH_Bohm string: ${base_input}$((${i}-1))${suffix}"
    echo "Replacement run_lammps_SPH_Bohm string: ${base_input}${i}${suffix}"
    sed -i -e "s/< ${base_input}$((${i}-1))${suffix}/< ${base_input}${i}${suffix}/g"  ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/ele_therm_research/run_lammps_SPH_Bohm
    
    echo "Start in.* string: ${start_deck_var}${i}${start_deck_let}"
    echo "Replacement in.* string: ${start_deck_var}${i}${suffix}"
    cp ${base_input}${i}a ${base_input}${i}${suffix}
    sed -i -e "s/${start_deck_var}${i}${start_deck_let}/${start_deck_var}${i}${suffix}/g" ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/ele_therm_research/${base_input}${i}${suffix}

    sbatch run_lammps_SPH_Bohm
    #cp ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/Al_research/3_5eV_5_2gcc/run_lammps_SPH_Bohm ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/Al_research/3_5eV_5_2gcc/run_lammps_SPH_Bohms${i}
    final_out=${base_input}${i}${suffix}
    fi
done    

sed -i -e "s/< ${final_out}/< in.placeholder/g" ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/ele_therm_research/run_lammps_SPH_Bohm