import csv
import os


def filter(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        for row in reader:
            if float(row[4]) < 0.05:
                try:
                    float(row[-2]) 
                    float(row[-1])  
                    species = points_rem(row[0])
                    if species.count('\t') == 6:
                        writer.writerow([species.replace('"', '')])  
                except ValueError:
                    continue


def points_rem(string):
    return string.replace('.', '\t')

def remove_quotes(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    content = content.replace('"', '')
    with open(file_path, 'w') as file:
        file.write(content)

def add_taxonomic_levels(file_path):
    taxonomic_levels = "Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n"
    with open(file_path, 'r') as file:
        content = file.read()
    content = taxonomic_levels + content
    with open(file_path, 'w') as file:
        file.write(content)
        
def run_file(input_file_name, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    output_tsv = os.path.join(output_folder, os.path.basename(input_file_name))
    filter(input_file_name, output_tsv)
    add_taxonomic_levels(output_tsv)
    remove_quotes(output_tsv)



input_folder = 'Output/LEFSE/step2v39'
output_folder = 'Output/DAS_ONLY_LEFSE/final_39'


for filename in os.listdir(input_folder):
    if filename.endswith('.res'):
        input_file_path = os.path.join(input_folder, filename)
        run_file(input_file_path, output_folder)

print("Elaboration complete!")
