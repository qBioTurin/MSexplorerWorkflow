#!/bin/bash

process_file() {
    local input_file="$1"
    local output_basename="$2"

    python3 /Scripts/create_dict.py "$input_file" "$folderHalf/${output_basename}.tsv"
    lefse-format_input.py "$input_file" "$folder1/${output_basename}.in" -c 1 -s -1 -u 2 -o 1000000
    run_lefse.py -l 2 "$folder1/${output_basename}.in" "$folder2/${output_basename}.res"
    python3 /Scripts/rem_line.py "$folder2/${output_basename}.res" "$folder3/${output_basename}.res"
    python3 /Scripts/fix_sp.py "$folderHalf/${output_basename}.tsv" "$folder3/${output_basename}.res" "$folder4/${output_basename}.res"

    echo "Elaboration of file $input_file complete"
}

input="$1"
path="/input_files/$input"
folderHalf="/dict"
folder1="step1_lefse"
folder2="input_files/unprocesssed_lefse"
folder3="lefse_raw_names"
folder4="/input_files/final_output_lefse"

mkdir -p "$folderHalf" "$folder1" "$folder2" "$folder3" "$folder4"

if [ -f "$path" ]; then
    BASENAME=$(basename "$input")
    process_file "$path" "$BASENAME"
elif [ -d "$path" ]; then
    for FILE in "$path"/*; do
        if [ -f "$FILE" ]; then
            echo "Elaboration of file: $FILE"
            BASENAME=$(basename "$FILE")
            process_file "$FILE" "$BASENAME"
        fi
    done
else
    echo "Error: $input is nor a file nor a directory."
    exit 1
fi
