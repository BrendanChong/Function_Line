#!/bin/bash

# Input file template
FILE="Single_Element.cpp"

# Get base name of the file (without extension)
BASENAME=$(basename "$FILE" .cpp)

rm -r *.txt
rm -r *.exo*
rm -r log_*

# source ~/.cshrc_moris

# Loop through restart files
for((R=1; R<=200; R+=1)) do
    # Use sed to replace the values
    # -i handles in-place editing
    sed -i "s/uint tNE = .*;/uint tNE = $R;/" "$FILE"

    echo "Generated: $FILE"

    moris.sh opt 1 "$BASENAME" 2>&1 | tee "log_${R}"

    # mv "IQI_Dof_Sensitivity.txt" "IQI_Dof_Sensitivity_${R}.txt"
    rm "IQI_Dof_Sensitivity.txt"
    rm "SEN.hdf5"
    mv "dv_dp.txt" "dv_dp_${R}.txt"
done
