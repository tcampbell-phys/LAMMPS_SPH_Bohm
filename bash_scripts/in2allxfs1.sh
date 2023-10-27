#!/bin/bash


rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../input_decks/allxfs1/Al_research/ campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/SPH_Bohm/Al_research/


# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxfs1/XC_tests/in.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/SPH_Bohm/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxfs1/small_tests/in.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/SPH_Bohm/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxfs1/dynamic_bohm/ion_electron_inputs/in.hedp_SPH_dynamic_bohm_21* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/SPH_Bohm/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxfs1/dynamic_bohm/ion_electron_inputs/in.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/SPH_Bohm/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxfs1/dynamic_bohm/bohm_only_inputs/in.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/SPH_Bohm/

#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../hedp_bohm_testing/in.pair_* campbell@allxfs1.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/pair_bohm_tests/
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/atom_vec_forceterms.* campbell@allxfs1.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/atom.* campbell@allxfs1.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  

#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/compute* campbell@allxfs1.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/dump_custom.* campbell@allxfs1.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/pair_bohm_global.* campbell@allxfs1.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/pair_bohm_global_ions.* campbell@allxfs1.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/atom_vec_SPH.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/atom.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/compute_* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/dump_custom.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_hybrid.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_cut_SPH.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_cut_pseudo.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_long_pseudo.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../../LAMMPS_Python/pseudopotentials/2LAMMPS/al.lda.lps.table campbell@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/SPH_Bohm/  
# scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../../LAMMPS_Python/pseudopotentials/2LAMMPS/al.lda.lps_even_temper.table campbell@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/SPH_Bohm/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_long_SPH.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_pseudo/pair_coul_long_SPH_ph.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_XC/* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_sub_XC/* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_bohm/* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_no_pseudo/* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_pseudo/* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_pseudo/pair_coul_cut_SPH.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/legacy_SPH_code/pair_coul_long_SPH_nopseudo.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/legacy_SPH_code/pair_coul_cut_SPH_nopseudo.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/update.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/KSPACE/pair_coul_long.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/KSPACE/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/KSPACE/ewald.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/KSPACE/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/pair.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/pair_coul_cut.cpp campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/compute_temp.cpp campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/fix_nve.cpp campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  


# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_bohm/pair_bohm_SPH_dynamic_Mocz* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/legacy_SPH_code/pair_bohm_SPH_dynamic_Mocz* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_dynamic_widths.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/min.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_fermi_no_cons.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_nh.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/thermo.* campbellt@allxfs1.physics.ox.ac.uk:~/fresh-lammps/lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/pair_bohm_global.* campbell@allxfs1.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  


# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ~/Documents/StructureFactorScripts/ campbell@allxfs1.physics.ox.ac.uk:~/DSF/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ~/Documents/StructureFactorScripts/* campbell@allxfs1.physics.ox.ac.uk:/data/tcampbell/StructureFactorScripts/

#ssh -J campbellt@bastion.physics.ox.ac.uk campbell@allxfs1.physics.ox.ac.uk 