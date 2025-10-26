# LAMMPS_SPH_Bohm

This repo contains the majority of Bohm SPH work conducted using a modified version of LAMMPS (version 3Mar20). Discussions of some of the results can be found at:
- Thesis: https://ora.ox.ac.uk/objects/uuid:4b9e073d-a4d3-47b7-b2cd-243bcc8d61a2/files/dsq87bv35z
- Paper: https://journals.aps.org/prresearch/abstract/10.1103/1gwd-z8zq
  

Folders:
- bash_scripts: bash code used to port data between various machines at Oxford (allxfs1, hedp64 (now allxdaq12)), Harwell (SCARF) and my local machine.
- input_decks: all LAMMPS input scripts used to run simulations on the various machines. The outputs of these scripts, minus the large dump files, can be found in the ./outputs folder. I rarely used allxfs1 so I doubt there is much interesting in here. hedp64 I used a lot for testing new implementations and performing sanity checks, as well as some small scale versions of the full simulations ran on SCARF. SCARF inputs contain the biggest simulations and is where I generated the majority of the published data.
- LAMMPS_code: the 3Mar20 LAMMPS src code and my additional bespoke folder USER-TC_SPH. Evertyhing required to build Bohm SPH within LAMMPS is in USER-TC_SPH.
- 


To start:
As stated, the Bohm SPH code is contained in the folder LAMMPS_code/src/USER-TC_SPH/. To get this to compile within LAMMPS, the individual files within USER-TC_SPH and it's subfolders need to be placed directly into the /src/ folder of LAMMPS on your target machine. Use the contents of atom_vec_SPH_CoM/ ahead of atom_vec_SPH/; atom_vec_SPH_CoM/ includes the variables required for the centre-of-mass formulation of Bohm SPH. Otherwise place all files (minus those contained in legacy_SPH_code/ and any duplicates of those in atom_vec_SPH_CoM/) into your src file and compile. Then, provided all the required libraries are available on the machine, compile with code such as: (from one of my compilation scripts):

"#!/bin/bash

echo "Cores: "
read cores

cd ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/build/

cmake ../cmake

cmake -D PKG_USER-MISC=yes -D PKG_KSPACE=on -D PKG_MANYBODY=yes -D LAMMPS_MACHINE=SPH_Bohm_191124_2 --build .

make -j ${cores} install"

Helpful to consult the manual on the inner workings of LAMMPS and tips on its compilation: https://docs.lammps.org/Manual.html

I recommend trying to re-run some of the input scripts used to generate outputs for the paper, these call the most up-to-date Bohm and SPH routines in the code base. In particular:
- input_decks/hedp64/CoM_thermo_tests/QHO_tests_241024/in.hedp_SPH_QHO_stability_test_20a to input_decks/hedp64/CoM_thermo_tests/QHO_tests_241024/in.hedp_SPH_QHO_stability_test_27a for Bohm conservation checks
- input_decks/hedp64/coul_conservation_tests/pre311024/in.hedp_SPH_coul_stability_test_2e_base and input_decks/hedp64/coul_conservation_tests/pre311024/in.hedp_SPH_coul_stability_test_2e_thetacut for the coulomb conservation checks
- 



I would ignore the Aluminium and Beryllium input scripts, there is an undiagnosed error in the energy conservation of the pseudopotential routine designed for the Gaussian SPH kernels.


This repo is not the easiest to navigate, for which I apologise in advance. In particular, there is no real separation of various versions of the code - I improved as I went. My ambitions of tidying the repo and writing an exhaustive guide remain unfulfilled, but I hope there is some value in making it all available. If hopelessly stuck you can email me at thomas_campbell1@protonmail.com and I will reply if I can.

Tom Campbell 26/10/25
