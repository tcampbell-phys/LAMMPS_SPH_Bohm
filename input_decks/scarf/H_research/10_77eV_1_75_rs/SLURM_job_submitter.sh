#!/bin/bash

### SLURM job submissions

start_string="in.placeholder"
base_input="in.scarf_H_10_77eV_1_75_rs_min_"

echo "Number range start: "
read start

echo "Number range end: "
read end

for ((i=$start; i<=$end; i++))
do
    if [ $i -eq ${start} ];
    then
    echo "Start string: ${start_string}"
    echo "Replacement string: ${base_input}${i}"
    sed -i -e "s/< ${start_string}/< ${base_input}${i}/g" ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/H_research/10_77eV_1_75_rs/run_lammps_SPH_Bohm
    sbatch run_lammps_SPH_Bohm
    final_out=${base_input}${i}
    #cp ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/H_research/10_77eV_1_75_rs/run_lammps_SPH_Bohm ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/H_research/10_77eV_1_75_rs/run_lammps_SPH_Bohm${i}
    fi
    if [ $i -ne ${start} ];
    then
    echo "Start string: ${base_input}$((${i}-1))"
    echo "Replacement string: ${base_input}${i}"
    sed -i -e "s/< ${base_input}$((${i}-1))/< ${base_input}${i}/g" ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/H_research/10_77eV_1_75_rs/run_lammps_SPH_Bohm
    sbatch run_lammps_SPH_Bohm
    #cp ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/H_research/10_77eV_1_75_rs/run_lammps_SPH_Bohm ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/H_research/10_77eV_1_75_rs/run_lammps_SPH_Bohms${i}
    final_out=${base_input}${i}
    fi
done

sed -i -e "s/< ${final_out}/< in.placeholder/g" ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/H_research/10_77eV_1_75_rs/run_lammps_SPH_Bohm

