#!/bin/bash

echo "Hedp number range start: "
read start
echo "Hedp number range end: "
read end
# echo "Letter variable (0 = a, 1 = b, 2 = c): "
# read letter_var

# if [ $letter_var -eq 0 ];
# then
#     suffix="a"
# fi

# if [ $letter_var -eq 1 ];
# then
#     suffix="b"
# fi

# if [ $letter_var -eq 2 ];
# then
#     suffix="c"
# fi

for ((i=$start; i<=$end; i++))
do

    #mkdir ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/hedp64/hedp_legacy_bohm_$i
    # mkdir ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/hedp64/static_width/hedp_SPH_bohm_$i
    # mkdir ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/hedp64/dynamic_width/hedp_SPH_dynamic_bohm_$i
    mkdir ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/hedp64/dynamic_width/hedp_SPH_small_test_$i


    # only use this command when target folder is small
    # scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/SPH_Bohm/hedp_legacy_bohm_$i${suffix}/*.dump ~/Mirror/LACIE_doc/LAMMPS/lammps-3Mar20/full_bohm_run/hedp_bohm_testing/hedp_bohm_testing_$i${suffix}/

    #scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/SPH_Bohm/hedp_legacy_bohm_$i/log.* ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/hedp64/hedp_legacy_bohm_$i/
    
    # scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/SPH_Bohm/hedp_SPH_bohm_$i/log.* ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/hedp64/static_width/hedp_SPH_bohm_$i/
    # scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/SPH_Bohm/hedp_SPH_bohm_$i/*.dump ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/hedp64/static_width/hedp_SPH_bohm_$i/
    
    # scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/SPH_Bohm/hedp_SPH_dynamic_bohm_$i/log.* ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/hedp64/dynamic_width/hedp_SPH_dynamic_bohm_$i/
    # scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/SPH_Bohm/hedp_SPH_dynamic_bohm_$i/*.dump ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/hedp64/dynamic_width/hedp_SPH_dynamic_bohm_$i/
    
    scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/SPH_Bohm/hedp_SPH_small_test_$i/log.* ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/hedp64/dynamic_width/hedp_SPH_small_test_$i/
    scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/SPH_Bohm/hedp_SPH_small_test_$i/*.dump ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/hedp64/dynamic_width/hedp_SPH_small_test_$i/
    
    # scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/SPH_Bohm/hedp_SPH_dynamic_bohm_$i/pair_coul_cut_SPH.* ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/hedp64/dynamic_width/hedp_SPH_dynamic_bohm_$i/

done