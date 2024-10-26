#!/bin/bash

### Systematic input script editor

# input coupling strengths
myArray=("0.2807897863750459" "0.5022513896051121" "0.7237134656383087" "0.945175143592007" "1.1666369912373342")
# base script template
template="scarf_e_therm_Z_4_coupling_scan_Nepe_8_"
letter="a"

count=12
for str in ${myArray[@]}; do
    # change base inputs to target values
    sed -i -e "s/string temp/string ${count}${letter}/g" /Users/tomcampbell/Documents/GitHub/LAMMPS_SPH_Bohm/input_decks/scarf/Electron_therm_research/Al_Z_4_coupling_scan/in.${template}temp
    sed -i -e "s/placeholder/$str/g" /Users/tomcampbell/Documents/GitHub/LAMMPS_SPH_Bohm/input_decks/scarf/Electron_therm_research/Al_Z_4_coupling_scan/in.${template}temp
    
    # copy to new file
    cp /Users/tomcampbell/Documents/GitHub/LAMMPS_SPH_Bohm/input_decks/scarf/Electron_therm_research/Al_Z_4_coupling_scan/in.${template}temp /Users/tomcampbell/Documents/GitHub/LAMMPS_SPH_Bohm/input_decks/scarf/Electron_therm_research/Al_Z_4_coupling_scan/in.${template}${count}${letter}
    
    echo "temp changed to ${count}${letter}"
    echo "Replacement string: $str"

    # return base script to template
    sed -i -e "s/string ${count}${letter}/string temp/g" /Users/tomcampbell/Documents/GitHub/LAMMPS_SPH_Bohm/input_decks/scarf/Electron_therm_research/Al_Z_4_coupling_scan/in.${template}temp
    sed -i -e "s/$str/placeholder/g" /Users/tomcampbell/Documents/GitHub/LAMMPS_SPH_Bohm/input_decks/scarf/Electron_therm_research/Al_Z_4_coupling_scan/in.${template}temp
    
    count=$((count+1))
done
