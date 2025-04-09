#!/bin/bash
folder1="Output/DAS_ONLY_LEFSE_HD/264_V1"
folder2="Output/DAS_ONLY_LEFSE_HD/264_V2"
DIRECTORY="Output/DAS_ONLY_LEFSE_HD/264_V0"

mkdir -p "$folder1"
mkdir -p "$folder2"

for FILE in "$DIRECTORY"/*; do  
    if [ -f "$FILE" ]; then
        echo "Processing file: $FILE"

        BASENAME=$(basename "$FILE")

        lefse-format_input.py "$FILE" "$folder1/${BASENAME}.in" -c 1 -s -1 -u 2 -o 1000000
        run_lefse2.py -w 1 -l 2 "$folder1/${BASENAME}.in" "$folder2/${BASENAME}.res"

        echo "Elaboration complete for $FILE"
    fi
done

