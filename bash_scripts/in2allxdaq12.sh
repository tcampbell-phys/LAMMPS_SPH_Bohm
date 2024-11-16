#!/bin/bash
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/multi_analyse_structure_fast_v2.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/test_sub_multi_analyse_structure_fast.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/structure_analysis_scripts/sub_multi_analyse_structure_fast.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/structure_analysis_scripts/ion_electron_structure.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/structure_analysis_scripts/CoM_check.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/thermalisation/therm_analysis.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/thermalisation/read_package/ campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/


# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/test_multi_analyse_structure_fast.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/analyse_structure_fast_v2.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/svensson_lib_v2/ campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/svensson_lib_v2/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/DSF_lib/ campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/DSF_lib/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/lammps_python_lib/ campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/lammps_python_lib/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/analyse_structure.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/analyse_structure_fast.py campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/svensson_lib/ campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/svensson_lib/
# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' ../../LAMMPS_Python/activate_python_env.sh campbell@allxdaq12.nat.physics.ox.ac.uk:~/structure_computation/

# rsync -azv -e 'ssh -o "ProxyCommand ssh -A campbellt@bastion.physics.ox.ac.uk -W %h:%p"' /Users/tomcampbell/Documents/GitHub/LAMMPS_SPH_Bohm/ref_outputs/AtomicKinetics/Al_research/YS_physical_ref.test campbell@allxdaq12.nat.physics.ox.ac.uk:/data/tcampbell/LAMMPS/SCARF_outputs/thermalisation/


# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxdaq12/XC_tests/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxdaq12/small_tests/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxdaq12/dynamic_bohm/ion_electron_inputs/in.hedp_SPH_dynamic_bohm_21* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxdaq12/dynamic_bohm/ion_electron_inputs/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/allxdaq12/dynamic_bohm/bohm_only_inputs/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/coul_long_sub_tests/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_coul_tests/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/confine_tests/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_confine/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/CoM_thermo_tests/QHO_ground_zeta_scan/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/CoM_thermo/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/CoM_thermo_tests/QHO_tests_241024/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/CoM_thermo/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/coul_conservation_tests/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/coul_conserve_tests/

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/coul_conservation_tests/081124/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/coul_conserve_tests/

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/coul_conservation_tests/pre311024/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/coul_conserve_tests/

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/coul_conservation_tests/coul_cut/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/coul_conserve_tests/

scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/coul_conservation_tests/coul_long/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/coul_conserve_tests/

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/coul_conservation_tests/small_system_tests/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/coul_conserve_tests/

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/MF_Bohm_tests/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/MF_Bohm_tests/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/H_ground_tests/current/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/H_ground_tests/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/BohmPotTests/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/BohmPotTests/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/Electron_therm_research/Al_Z_1_batch1/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/ele_therm_tests/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/Electron_therm_research/Al_Z_4_batch1/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/ele_therm_tests/

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/trajectories_tests/in.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/trajectories_tests/

#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../hedp_bohm_testing/in.pair_* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/pair_bohm_tests/
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/atom_vec_forceterms.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/atom.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  

#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/compute* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/dump_custom.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/pair_bohm_global.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/pair_bohm_global_ions.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/atom_vec_SPH.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/atom_vec_SPH_CoM/atom_vec_SPH_CoM.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/atom.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/compute_* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/dump_custom.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_hybrid.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/update.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_dynamic_widths.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_spherical_init.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/legacy_SPH_code/fix_dynamic_widths_CoM.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_dynamic_widths_CoM_MF.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_dynamic_widths_CoM_F2.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_dynamic_widths_CoM.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_dynamic_widths_CoM_noPBC.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_bohm_pot_calc_Mocz_robust_ftwo.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_confine/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_bohm_QHD/bohm_pot_in_fix/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_BP/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_dynamic_widths_CoM_BP_Mocz_rob.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  



# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_fermi_no_cons.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_nh_CoM.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_nvt_CoM.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  

### Trajectories Modules
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/KSPACE_traj/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_sub_XC_traj/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_bohm_QHD/pair_bohm_SPH_dynamic_Mocz_robust_F2_traj.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_sub_coul_no_pseudo/pair_coul_long_SPH_nopseudo_ph_sub_traj.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  


# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_cut_SPH.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_cut_pseudo.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_long_pseudo.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../../LAMMPS_Python/pseudopotentials/2LAMMPS/al.lda.lps.table campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/  
# scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../../LAMMPS_Python/pseudopotentials/2LAMMPS/al.lda.lps_even_temper.table campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_long_SPH.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_pseudo/pair_coul_long_SPH_ph.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_XC/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_confine/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_fixed_confine/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/legacy_SPH_code/pair_coul_long_SPH_nopseudo_ph_sub_leg.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_bohm_QHD/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_bohm_MF/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_no_pseudo/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_sub_coul_no_pseudo/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_sub_coul_no_pseudo/pair_coul_long_SPH_exp* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_sub_coul_no_pseudo_cut/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_sub_coul_pseudo/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_sub_coul_no_pseudo_MF/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_pseudo/* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_pseudo/pair_coul_cut_SPH.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/legacy_SPH_code/pair_coul_long_SPH_nopseudo.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/legacy_SPH_code/pair_coul_cut_SPH_nopseudo.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_sub_XC/pair_sub_Minoo* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/KSPACE/pair_coul_long.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/KSPACE/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/KSPACE/ewald.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/KSPACE/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/pair.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_nh.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_nvt_CoM.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/compute_temp_CoM.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/compute_temp.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/thermo.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_init_sample_vel.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_nve_BO.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/compute.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/pair_coul_cut.cpp campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/compute_temp.cpp campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/fix_nve.cpp campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  


# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_bohm/pair_bohm_SPH_dynamic_Mocz* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/legacy_SPH_code/pair_bohm_SPH_dynamic_Mocz* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/min.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_nh.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/thermo.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/pair_bohm_global.* campbell@allxdaq12.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  


# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ~/Documents/StructureFactorScripts/ campbell@allxdaq12.nat.physics.ox.ac.uk:~/DSF/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ~/Documents/StructureFactorScripts/* campbell@allxdaq12.nat.physics.ox.ac.uk:/data/tcampbell/StructureFactorScripts/

#ssh -J campbellt@bastion.physics.ox.ac.uk campbell@allxdaq12.nat.physics.ox.ac.uk 