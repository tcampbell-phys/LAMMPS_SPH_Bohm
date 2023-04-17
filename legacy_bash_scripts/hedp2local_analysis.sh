#!/bin/bash

echo "${prefix} folder (0 = ~/DSF/, 1= ~/DSF/StructureFactorScripts):"
read out

if [ $out -eq 0 ];

then
    folder="~/DSF/"
fi

if [ $out -eq 1 ];
then
    folder="/data/tcampbell/StructureFactorScripts/"
fi

echo "full_bohm run number: "
read number

echo "output (0 = Al, 1= Be):"
read out2

if [ $out2 -eq 0 ];
then
    outputdir="/Volumes/LaCie/Toolkit/Remote\ Files/StructureFactorFiles/StructureFactorScripts/"
    #mkdir /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/full_bohm_global_ion_al_${number}

    #outputdir="./Aluminium/"
    #mkdir ${outputdir}full_bohm_global_ion_al_${number}

    #scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:${folder}full_bohm_global_ion_al_${number}_data_Sk* ${outputdir}full_bohm_global_ion_al_${number}/
    scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:${folder}full_bohm_global_ion_al_${number}_electrons_* ${outputdir}
    #scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:${folder}StructureFactor_background.m /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/full_bohm_global_ion_be_${number}/
    #scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:${folder}DynamicStructureFactor_background.m /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/full_bohm_global_ion_be_${number}/

    #scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:${folder}full_bohm_global_ion_al_${number}_data_* /Volumes/LaCie/Toolkit/Remote\ Files/StructureFactorFiles/StructureFactorScripts/

fi

if [ $out2 -eq 1 ];
then

    outputdir="/Volumes/LaCie/Toolkit/Remote Files/StructureFactorFiles/StructureFactorScripts/"
    echo ${outputdir}
    #outputdir="/Volumes/LaCie/Toolkit/Remote Files/LAMMPS_backups/Beryllium/"
    #mkdir /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/full_bohm_global_ion_be_${number}

    #outputdir = './Beryllium/'
    #mkdir ${outputdir}full_bohm_global_ion_be_${number}

    scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:${folder}full_bohm_global_ion_be_${number}* ${outputdir}
    
    #scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:${folder}full_bohm_global_ion_be_${number}_data_* /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/full_bohm_global_ion_be_${number}/
    #scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:${folder}StructureFactor_background.m /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/full_bohm_global_ion_be_${number}/
    #scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:${folder}DynamicStructureFactor_background.m /Volumes/LaCie/Toolkit/Remote\ Files/LAMMPS_backups/full_bohm_global_ion_be_${number}/

    #scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" campbell@hedp64.nat.physics.ox.ac.uk:${folder}full_bohm_global_ion_be_${number}_data_* /Volumes/LaCie/Toolkit/Remote\ Files/StructureFactorFiles/StructureFactorScripts/

fi


#mkdir MC_run104_manual_tables_overlay
#scp -r campbell@hedp64.nat.physics.ox.ac.uk:~/DSF/MC_run104_manual_tables_overlay_data_* MC_run104_manual_tables_overlay/