#!/bin/bash

### SLURM job submissions

start_string="placeholder"
# input_files=("scarf_e_therm_Z_4_constant_KE_Nepe_scan_4a" "scarf_e_therm_Z_4_constant_KE_Nepe_scan_4b" "scarf_e_therm_Z_4_constant_KE_Nepe_scan_4c" "scarf_e_therm_Z_4_constant_KE_Nepe_scan_32a" "scarf_e_therm_Z_4_constant_KE_Nepe_scan_32b" "scarf_e_therm_Z_4_constant_KE_Nepe_scan_32c" "scarf_e_therm_Z_4_constant_KE_Nepe_scan_32a" "scarf_e_therm_Z_4_constant_KE_Nepe_scan_32b" "scarf_e_therm_Z_4_constant_KE_Nepe_scan_32c" "scarf_e_therm_Z_4_constant_KE_Nepe_scan_32a" "scarf_e_therm_Z_4_constant_KE_Nepe_scan_32b" "scarf_e_therm_Z_4_constant_KE_Nepe_scan_32c")

input_files=("scarf_e_therm_Z_4_constant_KE_Nepe_scan_32a" "scarf_e_therm_Z_4_constant_KE_Nepe_scan_32b" "scarf_e_therm_Z_4_constant_KE_Nepe_scan_32c" "scarf_e_therm_Z_4_constant_temp_Nepe_scan_32a" "scarf_e_therm_Z_4_constant_temp_Nepe_scan_32b" "scarf_e_therm_Z_4_constant_temp_Nepe_scan_32c")

for str in ${input_files[@]}; do
  echo "Start run_lammps_SPH_Bohm string: ${start_string}"
  echo "Replacement run_lammps_SPH_Bohm string: $str"
  sed -i -e "s/${start_string}/${str}/g" ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/ele_therm_research/run_lammps_SPH_Bohm
  sbatch run_lammps_SPH_Bohm
  sed -i -e "s/${str}/${start_string}/g" ~/mywork/LAMMPS/fresh-lammps-3Mar20/lammps-3Mar20/ele_therm_research/run_lammps_SPH_Bohm
done