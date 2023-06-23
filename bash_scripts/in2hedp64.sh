#!/bin/bash

scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../input_decks/hedp64/dynamic_bohm/in.* campbell@hedp64.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/SPH_Bohm/
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../hedp_bohm_testing/in.pair_* campbell@hedp64.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/pair_bohm_tests/
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/atom_vec_forceterms.* campbell@hedp64.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/atom.* campbell@hedp64.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  

#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/compute* campbell@hedp64.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/dump_custom.* campbell@hedp64.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/pair_bohm_global.* campbell@hedp64.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/pair_bohm_global_ions.* campbell@hedp64.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/* campbell@hedp64.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/atom_vec_SPH.* campbell@hedp64.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/atom.* campbell@hedp64.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/compute_* campbell@hedp64.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/dump_custom.* campbell@hedp64.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_coul_cut_SPH.* campbell@hedp64.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/pair_bohm_SPH_dynamic_Mocz.* campbell@hedp64.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/legacy_SPH_code/pair_bohm_SPH_dynamic_Mocz* campbell@hedp64.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/fix_dynamic_widths.* campbell@hedp64.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/src/USER-TC_SPH/thermo.* campbell@hedp64.nat.physics.ox.ac.uk:~/fresh_lammps-3Mar20/lammps-3Mar20/src/  

# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ../LAMMPS_code/pair_bohm_global.* campbell@hedp64.nat.physics.ox.ac.uk:~/allxfs1-lammps-3Mar20/src/  


# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ~/Documents/StructureFactorScripts/ campbell@hedp64.nat.physics.ox.ac.uk:~/DSF/
# scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ~/Documents/StructureFactorScripts/* campbell@hedp64.nat.physics.ox.ac.uk:/data/tcampbell/StructureFactorScripts/

