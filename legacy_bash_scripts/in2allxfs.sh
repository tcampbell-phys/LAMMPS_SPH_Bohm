#!/bin/bash

#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" ./Beryllium/in.* campbellt@allxfs1.physics.ox.ac.uk:~/lammps-3Mar20/pair_bohm_global/
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" pauli*.table campbellt@allxfs1.physics.ox.ac.uk:~/lammps-3Mar20/pair_bohm_global/
scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" in.pauli_scan_6 campbellt@allxfs1.physics.ox.ac.uk:~/lammps-3Mar20/pair_bohm_global/

#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" Yukawa_Tests/in.aluminum* campbellt@allxfs1.physics.ox.ac.uk:~/lammps-3Mar20/Yukawa/
#scp -r -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" pseudo*.table campbellt@allxfs1.physics.ox.ac.uk:~/lammps-3Mar20/pair_bohm_global/
#scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" zero.table campbellt@allxfs1.physics.ox.ac.uk:~/lammps-3Mar20/pair_bohm_global/
#scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" full_pauli.table campbellt@allxfs1.physics.ox.ac.uk:~/lammps-3Mar20/pair_bohm_global/

#scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" pair_bohm_global*.cpp campbellt@allxfs1.physics.ox.ac.uk:~/lammps-3Mar20/src/
#scp -oProxyCommand="ssh -W %h:%p campbellt@bastion.physics.ox.ac.uk" pair_bohm_global*.h campbellt@allxfs1.physics.ox.ac.uk:~/lammps-3Mar20/src/
