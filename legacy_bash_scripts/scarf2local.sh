#!/bin/bash
echo "${prefix} suffix (0 = 'scarf_pauli_scan_',1 = 'full_bohm_global_ion_be_',2 = 'full_bohm_global_ion_al_):"
read out
suffix=""

if [ $out -eq 0 ];
then
    folder="Bohm_Global"
    prefix="scarf_pauli_scan_"
    out_prefix=""
fi
if [ $out -eq 1 ];
then
    folder="Bohm_Global"
    prefix="full_bohm_global_ion_be_"
    out_prefix='Beryllium/'
fi
if [ $out -eq 2 ];
then
    folder="Bohm_Global"
    prefix="full_bohm_global_ion_al_"
    out_prefix='Aluminium/'
fi

echo "Number range start: "
read start

echo "Number range end: "
read end

echo "Letter variable (0 = a, 1 = b, 2 = c, 3 = d, 4 = e, 5 = none): "
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

if [ $letter_var -eq 4 ];
then
    suffix="e"
fi

if [ $letter_var -eq 5 ];
then
    suffix=""
fi

for ((i=$start; i<=$end; i++))
do

    mkdir ${out_prefix}${prefix}$i${suffix}

    scp scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/${folder}/${prefix}$i${suffix}/log.${prefix}$i${suffix} ./${out_prefix}${prefix}$i${suffix}
    # only have the below scp command on for very small all_data dump files
    scp scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/${folder}/${prefix}$i${suffix}/all_data.dump ./${out_prefix}${prefix}$i${suffix}

    #scp -r scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/${folder}/${prefix}$i${suffix}/*.txt ./${out_prefix}${prefix}$i${suffix}
    #scp -r scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/${folder}/${prefix}$i${suffix}/*.rdf ./${out_prefix}${prefix}$i${suffix}

    #scp -r scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/${prefix}$i/*.rdf ./${prefix}$i
    #scp -r scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/${folder}/${prefix}$i${suffix}/*txt* ./${prefix}_hedp$i${suffix}/
    
done