#!/bin/bash

#scp  data2hedp64.sh scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/Bohm_Global/bash_scripts/

#rsync -auv  ./Beryllium/in.pair_bohm_global_ion* scarf734@ui1.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/
#rsync -auv ./Aluminium/in.pair_bohm_global_ion* scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/lammps-3Mar20/pair_bohm_global/historic_inputs/
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
scp ../LAMMPS_code/src/USER-TC_SPH/*.cpp scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
scp ../LAMMPS_code/src/USER-TC_SPH/*.h scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# scp ../LAMMPS_code/src/USER-TC_SPH/pair_coul_pseudo/* scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# scp ../LAMMPS_code/src/USER-TC_SPH/pair_coul_no_pseudo/* scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/
# scp ../LAMMPS_code/src/USER-TC_SPH/pair_bohm/* scarf734@ui2.scarf.rl.ac.uk:~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/src/