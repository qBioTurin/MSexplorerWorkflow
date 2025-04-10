import csv
import re
import sys
"""Normalize name by removing non-alphabetic characters and converting to lowercase."""
def clean_name(name):
    return re.sub(r'[^a-zA-Z]', '', name).lower()

"""Create a dictionary with normalized names as keys and original names as values."""
def load_names(file_path):
    name_dict = {}
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            original_name = line.strip()
            cleaned_name = clean_name(original_name)
            name_dict[cleaned_name] = original_name
    return name_dict

"""Switch wrong name with the correct one ."""
def fix_taxonomy_file(taxonomy_file, correct_names, output_file):
    with open(taxonomy_file, 'r', encoding='utf-8') as f, open(output_file, 'w', encoding='utf-8', newline='') as out_f:
        reader = csv.reader(f, delimiter='\t')
        writer = csv.writer(out_f, delimiter='\t')
        
        header = next(reader)
        writer.writerow(header)
        
        for row in reader:
            species_name = row[-1] 
            cleaned_species = clean_name(species_name) 
            if cleaned_species in correct_names:
                row[-1] = correct_names[cleaned_species] 
            writer.writerow(row)


def main():
    if len(sys.argv) != 4:
        sys.exit(1)

    correct_names_file = sys.argv[1]
    taxonomy_file = sys.argv[2]
    output_file = sys.argv[3]

    correct_names = load_names(correct_names_file)
    fix_taxonomy_file(taxonomy_file, correct_names, output_file)
    print(f"Name correction complete")


if __name__ == "__main__":
    main()

