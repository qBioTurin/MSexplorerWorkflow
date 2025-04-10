import csv
import sys
import os

def filter(input_file, output_file):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter='\t')

            for row in reader:
                if len(row) >= 5:
                    try:
                        float(row[-2])  
                        float(row[-1])  
                        species = points_rem(row[0])
                        if species.count('\t') == 6:
                            writer.writerow([species.replace('"', '')])  
                    except ValueError:
                        continue
                else:
                    continue
    except Exception as e:
        print(f"Error during file processing: {e}")

def points_rem(string):
    transformed = string.replace('.', '\t')
    return transformed

def remove_quotes(file_path):
    try:
        with open(file_path, 'r') as file:
            content = file.read()
        content = content.replace('"', '')
        with open(file_path, 'w') as file:
            file.write(content)
    except Exception as e:
        print(f"Error removing quotes: {e}")

def add_taxonomic_levels(file_path):
    try:
        taxonomic_levels = "Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n"
        with open(file_path, 'r') as file:
            content = file.read()
        content = taxonomic_levels + content
        with open(file_path, 'w') as file:
            file.write(content)
    except Exception as e:
        print(f"Error adding taxonomic levels: {e}")

def run_file(input_file_name, output_file):
    try:
        output_tsv = output_file
        filter(input_file_name, output_tsv)
        add_taxonomic_levels(output_tsv)
        remove_quotes(output_tsv)
    except Exception as e:
        print(f"Error during run_file execution: {e}")

if len(sys.argv) != 3:
    print("Usage: python rem_line.py <input_file> <output_file>")
    sys.exit(1)

input_file_path = sys.argv[1]
output_file = sys.argv[2]
run_file(input_file_path, output_file)

print("Elaboration complete!")


