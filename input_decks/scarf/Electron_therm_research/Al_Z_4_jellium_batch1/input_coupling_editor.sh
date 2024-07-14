#!/bin/bash

### Systematic input script editor

# input coupling strengths
myArray=("1.2916459842148007" "1.4577423889591685" "1.6238387937035366" "1.789935198447905" "1.956031603192273" "2.122128007936641" "2.2882244126810094" "2.454320817425377" "2.6204172221697455" "2.7865136269141133" "2.952610031658482")

# base script template
template="scarf_e_therm_Z_4_gamma_scan_Nepe_16_jellium_"

count=1
for str in ${myArray[@]}; do
    # change base inputs to target values
    sed -i -e "s/string temp/string ${count}a/g" /Users/tomcampbell/Documents/GitHub/LAMMPS_SPH_Bohm/input_decks/scarf/Electron_therm_research/Al_Z_4_jellium_batch1/in.${template}temp
    sed -i -e "s/placeholder/$str/g" /Users/tomcampbell/Documents/GitHub/LAMMPS_SPH_Bohm/input_decks/scarf/Electron_therm_research/Al_Z_4_jellium_batch1/in.${template}temp
    
    # copy to new file
    cp /Users/tomcampbell/Documents/GitHub/LAMMPS_SPH_Bohm/input_decks/scarf/Electron_therm_research/Al_Z_4_jellium_batch1/in.${template}temp /Users/tomcampbell/Documents/GitHub/LAMMPS_SPH_Bohm/input_decks/scarf/Electron_therm_research/Al_Z_4_jellium_batch1/in.${template}${count}a
    
    echo "temp changed to ${count}a"
    echo "Replacement string: $str"

    # return base script to template
    sed -i -e "s/string ${count}a/string temp/g" /Users/tomcampbell/Documents/GitHub/LAMMPS_SPH_Bohm/input_decks/scarf/Electron_therm_research/Al_Z_4_jellium_batch1/in.${template}temp
    sed -i -e "s/$str/placeholder/g" /Users/tomcampbell/Documents/GitHub/LAMMPS_SPH_Bohm/input_decks/scarf/Electron_therm_research/Al_Z_4_jellium_batch1/in.${template}temp
    
    count=$((count+1))
done
