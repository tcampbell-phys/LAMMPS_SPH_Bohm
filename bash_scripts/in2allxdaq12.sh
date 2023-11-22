#!/bin/bash
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/multi_analyse_structure_fast_v2.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/analyse_structure_fast_v2.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/svensson_lib_v2/ campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/svensson_lib_v2/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/DSF_lib/ campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/DSF_lib/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/lammps_python_lib/ campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/lammps_python_lib/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/analyse_structure.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/analyse_structure_fast.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/svensson_lib/ campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/svensson_lib/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/activate_python_env.sh campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/



# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxdaq12/XC_tests/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxdaq12/small_tests/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxdaq12/dynamic_bohm/ion_electron_inputs/in.hedp_SPH_dynamic_bohm_21* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxdaq12/dynamic_bohm/ion_electron_inputs/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxdaq12/dynamic_bohm/bohm_only_inputs/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/

scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/coul_long_sub_tests/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_coul_tests/

#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../hedp_bohm_testing/in.pair_* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/pair_bohm_tests/
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/atom_vec_forceterms.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/atom.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  

#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/compute* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/dump_custom.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/pair_bohm_global.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/pair_bohm_global_ions.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/atom_vec_SPH.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/atom.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/compute_* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/dump_custom.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_cut_SPH.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_cut_pseudo.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_long_pseudo.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../../LAMMPS_Python/pseudopotentials/2LAMMPS/al.lda.lps.table campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/  
# scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../../LAMMPS_Python/pseudopotentials/2LAMMPS/al.lda.lps_even_temper.table campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_long_SPH.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_pseudo/pair_coul_long_SPH_ph.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_XC/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_bohm/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_no_pseudo/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_pseudo/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_pseudo/pair_coul_cut_SPH.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/legacy_SPH_code/pair_coul_long_SPH_nopseudo.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/legacy_SPH_code/pair_coul_cut_SPH_nopseudo.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/update.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/KSPACE/pair_coul_long.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/KSPACE/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/KSPACE/ewald.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/KSPACE/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/pair.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/pair_coul_cut.cpp campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/compute_temp.cpp campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/fix_nve.cpp campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  


# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_bohm/pair_bohm_SPH_dynamic_Mocz* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/legacy_SPH_code/pair_bohm_SPH_dynamic_Mocz* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_dynamic_widths.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/min.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_fermi_no_cons.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_nh.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/thermo.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/pair_bohm_global.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  


# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ~/Documents/StructureFactorScripts/ campbell@allxdaq12.nat.physics.ox.ac.uk:~/DSF/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ~/Documents/StructureFactorScripts/* campbell@allxdaq12.nat.physics.ox.ac.uk:/data/tcampbell/StructureFactorScripts/

#ssh -J campbellt@bastion.physics.ox.ac.uk campbell@allxdaq12.nat.physics.ox.ac.uk 