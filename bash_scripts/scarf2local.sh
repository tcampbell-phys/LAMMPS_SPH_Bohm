#!/bin/bash

echo "SCARF number range start: "
read start
echo "SCARF number range end: "
read end

for ((i=$start; i<=$end; i++))
do
    mkdir ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/SCARF/testing/scarf_SPH_full_test_$i


    # only use this command when target folder is small
    
    scp -r scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/SPH_Bohm/testing/scarf_SPH_full_test_$i/log.* ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/SCARF/testing/scarf_SPH_full_test_$i
    # scp -r scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/SPH_Bohm/testing/scarf_SPH_full_test_$i/*.dump ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/SCARF/testing/scarf_SPH_full_test_$i
    # scp -r scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/SPH_Bohm/testing/scarf_SPH_full_test_$i/*.rdf ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/SCARF/testing/scarf_SPH_full_test_$i

done