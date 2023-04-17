#!/bin/bash
echo "${prefix} suffix (0 = '_hedp_global_scan_', 1='_global_ion_al_', 2='pauli_scan_',3='_global_ion_',4='_global_ion_be_'):"
read out
suffix=""
if [ $out -eq 0 ];
then
    folder="Bohm_Global"
    prefix="full_bohm"
    out_prefix=""
    suffix="_hedp_global_scan_"
fi

if [ $out -eq 1 ];
then
    folder="Bohm_Global"
    prefix="full_bohm"
    out_prefix=""
    suffix="_global_ion_al_"
fi

if [ $out -eq 2 ];
then
    folder="Pauli_Scan"
    out_prefix=""
    prefix=""
    suffix="pauli_scan_"
fi

if [ $out -eq 3 ];
then
    folder="Bohm_Global"
    prefix="full_bohm"
    out_prefix=""
    suffix="_global_ion_"
fi

if [ $out -eq 4 ];
then
    folder="Bohm_Global"
    out_prefix="Beryllium/"
    prefix="full_bohm"
    suffix="_global_ion_be_"
    echo $out_prefix
fi

echo "run number: "
read number

if [ $out -eq 0 ];
then
    mkdir ${prefix}${suffix}${number}
fi

if [ $out -eq 1 ];
then
    mkdir ${prefix}${suffix}${number}
fi

if [ $out -eq 2 ];
then
    mkdir ${suffix}${number}
fi

if [ $out -eq 3 ];
then
    mkdir ${prefix}${suffix}${number}
fi

if [ $out -eq 4 ];
then
    mkdir "/Volumes/LaCie/Toolkit/Remote Files/LAMMPS_backups/"${out_prefix}${prefix}${suffix}${number}
fi

echo "output directory = ${out_prefix}${prefix}${suffix}${number}"

scp  -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbellt@allxfs1.physics.ox.ac.uk:/data/gregori/campbellt/LAMMPS/${folder}/${prefix}${suffix}${number}/log.${prefix}${suffix}${number} "/Volumes/LaCie/Toolkit/Remote Files/LAMMPS_backups/"${out_prefix}${prefix}${suffix}${number}
scp -r  -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbellt@allxfs1.physics.ox.ac.uk:/data/gregori/campbellt/LAMMPS/${folder}/${prefix}${suffix}${number}/*.txt "/Volumes/LaCie/Toolkit/Remote Files/LAMMPS_backups/"${out_prefix}${prefix}${suffix}${number}
scp -r  -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbellt@allxfs1.physics.ox.ac.uk:/data/gregori/campbellt/LAMMPS/${folder}/${prefix}${suffix}${number}/*.rdf "/Volumes/LaCie/Toolkit/Remote Files/LAMMPS_backups/"${out_prefix}${prefix}${suffix}${number}

#scp -r 1234 campbellt@127.0.0.1:/data/gregori/campbellt/LAMMPS//${prefix}${number}/*.rdf ./${prefix}${number}

#scp -r campbellt@allxfs1.physics.ox.ac.uk:/data/gregori/campbellt/LAMMPS/${prefix}${number}/*.rdf ./${prefix}${number}

#scp -r campbellt@allxfs1.physics.ox.ac.uk:/data/gregori/campbellt/LAMMPS/${folder}/${prefix}${suffix}${number}/*txt* ./${prefix}_hedp${suffix}${number}/
scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbellt@allxfs1.physics.ox.ac.uk:/data/gregori/campbellt/LAMMPS/${folder}/${prefix}${suffix}${number}/all_*.dump "/Volumes/LaCie/Toolkit/Remote Files/LAMMPS_backups/"${out_prefix}${prefix}${suffix}${number}
#scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbellt@allxfs1.physics.ox.ac.uk:/data/gregori/campbellt/LAMMPS/${folder}/${prefix}${suffix}${number}/all_*.dump "/Volumes/LaCie/Toolkit/Remote Files/LAMMPS_backups/"${prefix}${suffix}${number}
