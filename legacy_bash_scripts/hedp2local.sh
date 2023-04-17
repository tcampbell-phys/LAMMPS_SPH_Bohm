#!/bin/bash

echo "Hedp number range start: "
read start
echo "Hedp number range end: "
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

for ((i=$start; i<=$end; i++))
do

    mkdir ~/Mirror/LACIE_doc/LAMMPS/lammps-3Mar20/full_bohm_run/hedp_bohm_testing/hedp_bohm_testing_$i${suffix}

    # only use this command when target folder is small
    scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/Bohm_Global/hedp_bohm_testing_$i${suffix}/*.dump ~/Mirror/LACIE_doc/LAMMPS/lammps-3Mar20/full_bohm_run/hedp_bohm_testing/hedp_bohm_testing_$i${suffix}/

    scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/Bohm_Global/hedp_bohm_testing_$i${suffix}/log.* ~/Mirror/LACIE_doc/LAMMPS/lammps-3Mar20/full_bohm_run/hedp_bohm_testing/hedp_bohm_testing_$i${suffix}/
done