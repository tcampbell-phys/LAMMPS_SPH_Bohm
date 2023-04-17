#!/bin/bash
echo "${prefix} suffix (0 = 'scarf_pauli_scan_',1 = 'full_bohm_global_ion_be_',2 = 'full_bohm_global_ion_al_):"
read out
suffix=""

echo "${prefix} forward onto hedp64 (1 = yes):"
read forward
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

#for ((i=$start; i<=$end; i++))
#do
#    mkdir /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/${out_prefix}${prefix}$i${suffix}

#    scp scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/${folder}/${prefix}$i${suffix}/log.${prefix}$i${suffix} /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/${out_prefix}${prefix}$i${suffix}/
    #scp -r scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/${folder}/${prefix}$i${suffix}/*.txt /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/${out_prefix}${prefix}$i${suffix}/
    #scp -r scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/${folder}/${prefix}$i${suffix}/*.rdf /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/${out_prefix}${prefix}$i${suffix}/
#    scp scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/${folder}/${prefix}$i${suffix}/ions_xyz_data.dump /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/${out_prefix}${prefix}$i${suffix}/
#    scp -r scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/${folder}/${prefix}$i${suffix}/electrons_xyz_data* /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/${out_prefix}${prefix}$i${suffix}/
#done

if [ $forward -eq 1 ];
then
for ((i=$start; i<=$end; i++))
do
    scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/${out_prefix}${prefix}$i${suffix}/ions_xyz_data* campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/Bohm_Ions/${out_prefix}${prefix}$i${suffix}/
    scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/${out_prefix}${prefix}$i${suffix}/electrons_xyz_data* campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/Bohm_Ions/${out_prefix}${prefix}$i${suffix}/
done
fi
