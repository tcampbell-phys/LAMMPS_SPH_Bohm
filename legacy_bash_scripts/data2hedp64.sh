#!/bin/bash

echo "${prefix} prefix (0 = 'scarf_pauli_scan_',1 = 'full_bohm_global_ion_be_',2 = 'full_bohm_global_ion_al_):"
read out
suffix=""

if [ $out -eq 0 ];
then
    folder="Bohm_Global/Pauli_Scan"
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

echo "Letter variable (0 = a, 1 = b, 2 = c): "
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

#ssh -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/Bohm_Ions    
    
#for ((i=$start; i<=$end; i++))
#do
#    mkdir ${out_prefix}${prefix}$i${suffix}
#done

#~.

for ((i=$start; i<=$end; i++))
do
    cd ../${prefix}$i${suffix}
    scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" log.${prefix}$i${suffix}  campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/Bohm_Ions/${out_prefix}${prefix}$i${suffix}
    scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ions_xyz_data.dump  campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/Bohm_Ions/${out_prefix}${prefix}$i${suffix}
done