import csv
import sys
def extract_species(input_file, output_file):
    with open(input_file, 'r', encoding='utf-8') as infile, open(output_file, 'w', encoding='utf-8', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        next(reader)
        next(reader)
        
        species_names = set()  
        
        for row in reader:
            if row:  
                taxonomy = row[0]  
                species_name = taxonomy.split('|')[-1]  
                species_names.add(species_name)
        
        
        for species in sorted(species_names):  
            writer.writerow([species])


input_file = sys.argv[1]
output_file = sys.argv[2]
extract_species(input_file, output_file)
