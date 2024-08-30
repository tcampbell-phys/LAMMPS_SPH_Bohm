#!/bin/bash

### Systematic input script editor

# input coupling strengths
myArray=("1.1666370288543766" "1.3880989018468675" "1.609560774839358" "1.8310226478318492" "2.0524845208243403" "2.273946393816831" "2.4954082668093216" "2.7168701398018125" "2.9383320127943033" "3.159793885786794" "3.3812557587792846")
# base script template
template="scarf_e_therm_Z_4_coupling_scan_Nepe_8_"
letter="c"

count=1
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
