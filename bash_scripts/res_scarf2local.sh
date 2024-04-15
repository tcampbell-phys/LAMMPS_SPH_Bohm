#!/bin/bash

# rsync -av --exclude='*.dump' --exclude='*.cpp' --exclude='*.h' scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/SPH_Bohm/H_research/ ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/SCARF/H_research/
# rsync -av --exclude='*.dump' --exclude='*.cpp' --exclude='*.h' scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/SPH_Bohm/H_research/21_54eV_1_75_rs/ ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/SCARF/H_research/21_54eV_1_75_rs/

# rsync -av --exclude='*.dump' --exclude='*.cpp' --exclude='*.h' scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/SPH_Bohm/Al_research/ ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/SCARF/Al_research/
# rsync -av --exclude='*.dump' --exclude='*.cpp' --exclude='*.h' scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/SPH_Bohm/Be_research/ ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/SCARF/Be_research/

rsync -av --exclude='*.cpp' --exclude='*.h' scarf734@ui1.scarf.rl.ac.uk:/work4/clf/scarf734/LAMMPS/SPH_Bohm/ground/ ~/Documents/GitHub/LAMMPS_SPH_Bohm/outputs/SCARF/H_research/H_ground/
