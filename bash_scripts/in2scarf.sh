#!/bin/bash

#scp  data2hedp64.sh scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/Bohm_Global/bash_scripts/

#rsync -auv  ./Beryllium/in.pair_bohm_global_ion* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/
# rsync -auv ./Aluminium/in.pair_bohm_global_ion* scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/historic_inputs/
#rsync -auv  ./Pauli_Scan_Tests/in.scarf_pauli_* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/
#rsync -auv ../distribution_reader.py scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/
#scp -r  Yukawa_Tests/in.aluminum* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/Yukawa/
#rsync -auv  pauli*.table scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/
#rsync -auv  SLURM_job_submitter.sh scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/
#rsync -auv  pseudo*.table scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/
#scp  zero.table scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/
#scp  full_pauli.table scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/
#scp  fix_fermi_no_cons.cpp scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/src/
#scp  fix_fermi_no_cons.h scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/src/

# rsync ../input_decks/scarf/SPH_tests/in.* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/SPH_Bohm/
# rsync ../input_decks/scarf/Al_research/3_5eV_5_2gcc/* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/Al_research/3_5eV_5_2gcc/
# rsync ../input_decks/scarf/Al_research/0_5eV_2_7gcc/in.* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/Al_research/0_5eV_2_7gcc/
# rsync ../input_decks/scarf/Al_research/0_088eV_2_35gcc/in.* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/Al_research/0_088eV_2_35gcc/
# rsync ../input_decks/scarf/Al_research/1_1eV_3_4gcc/in.* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/Al_research/1_1eV_3_4gcc/
# rsync ../input_decks/scarf/Be_research/in.* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/Be_research/
rsync ../input_decks/scarf/H_research/15_71eV_1_26_rs/* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/H_research/15_71eV_1_26_rs/
rsync ../input_decks/scarf/H_research/12_53eV_2_rs/* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/H_research/12_53eV_2_rs/
rsync ../input_decks/scarf/H_research/10_77eV_1_75_rs/* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/H_research/10_77eV_1_75_rs/

# scp ../LAMMPS_code/src/USER-TC_SPH/*.cpp scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# scp ../LAMMPS_code/src/USER-TC_SPH/*.h scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# scp ../LAMMPS_code/src/USER-TC_SPH/pair_coul_pseudo/* scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# rsync ../LAMMPS_code/src/USER-TC_SPH/pair_coul_pseudo/* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# rsync ../LAMMPS_code/src/USER-TC_SPH/pair_sub_XC/* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# rsync ../LAMMPS_code/src/USER-TC_SPH/pair_coul_no_pseudo/* scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# scp ../LAMMPS_code/src/USER-TC_SPH/pair_bohm/* scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# scp ../LAMMPS_code/src/USER-TC_SPH/pair_XC/* scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# scp ../LAMMPS_code/src/USER-TC_SPH/update* scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# scp ../LAMMPS_code/src/USER-TC_SPH/fix_fermi_no_cons.* scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# scp ../LAMMPS_code/src/USER-TC_SPH/fix_dynamic_widths.* scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# scp ../LAMMPS_code/src/pair_coul_cut.cpp scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# scp ../LAMMPS_code/src/compute_temp.cpp scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# scp ../LAMMPS_code/src/fix_nve.cpp scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/