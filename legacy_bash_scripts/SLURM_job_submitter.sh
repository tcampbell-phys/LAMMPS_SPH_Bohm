#!/bin/bash

### SLURM job submissions

start_string="in.placeholder"
base_input="in.pair_bohm_global_ion_be_"

echo "Number range start: "
read start

echo "Number range end: "
read end

echo "Letter variable (0 = a, 1 = b, 2 = c): "
read letter_var

if [ $letter_var -eq 0 ];
then
    letter="a"
fi

if [ $letter_var -eq 1 ];
then
    letter="b"
fi

if [ $letter_var -eq 2 ];
then
    letter="c"
fi

echo ${letter}

for ((i=$start; i<=$end; i++))
do
    if [ $i -eq ${start} ];
    then
    echo "Start string: ${start_string}"
    echo "Replacement string: ${base_input}${i}${letter}"
    sed -i -e "s/< ${start_string}/< ${base_input}${i}${letter}/g" ~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/run_lammps_pair_bohm_global
    sbatch run_lammps_pair_bohm_global
    final_out=${base_input}${i}${letter}
    #cp ~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/run_lammps_pair_bohm_global ~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/run_lammps_pair_bohm_global${i}
    fi
    if [ $i -ne ${start} ];
    then
    echo "Start string: ${base_input}$((${i}-1))${letter}"
    echo "Replacement string: ${base_input}${i}${letter}"
    sed -i -e "s/< ${base_input}$((${i}-1))${letter}/< ${base_input}${i}${letter}/g" ~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/run_lammps_pair_bohm_global
    sbatch run_lammps_pair_bohm_global
    #cp ~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/run_lammps_pair_bohm_global ~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/run_lammps_pair_bohm_globals${i}
    final_out=${base_input}${i}${letter}
    fi
done

sed -i -e "s/< ${final_out}/< in.placeholder/g" ~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/run_lammps_pair_bohm_global

