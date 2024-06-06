import csv

# Define file paths
INPUT_FILE_STY = 'data/input/umls/2024AA/META/MRSTY.RRF'
INPUT_FILE_DEF = 'data/input/umls/2024AA/META/MRDEF.RRF'
INPUT_FILE_CONSO = 'data/input/umls/2024AA/META/MRCONSO.RRF'
OUTPUT_FILE_DISORDER = 'data/output/umls/umls_def_disorder.tsv'
OUTPUT_FILE_DISEASE = 'data/output/umls/umls_def_disease.tsv'
OUTPUT_FILE_CONSO = 'data/output/umls/umls.tsv'

# Define columns for MRSTY.RRF based on MRSTY.ctl
COLUMNS_STY = ['CUI', 'TUI', 'STN', 'STY', 'ATUI', 'CVF']

# Disorder entities
DISORDER_ENTITIES = [
    'Congenital Abnormality', 
    'Acquired Abnormality',
    'Injury or Poisoning', 
    'Pathologic Function', 
    'Disease or Syndrome', 
    'Mental or Behavioral Dysfunction', 
    'Cell or Molecular Dysfunction', 
    'Experimental Model of Disease', 
    'Sign and Symptom', 
    'Anatomical Abnormality', 
    'Neoplastic Process'
]

# English sources for MRDEF.RRF
ENGLISH_SOURCES = [
    'MSH', 'CSP', 'NCI', 'PDQ', 'NCI_NCI-GLOSS', 'CHV', 'NCI_CRCH', 'NCI_CareLex', 'NCI_CDISC-GLOSS', 
    'UWDA', 'FMA', 'SNOMEDCT_US', 'HPO', 'NCI_CTCAE', 'NCI_NICHD', 'MEDLINEPLUS', 'NCI_ACC-AHA', 
    'NCI_CDISC', 'NCI_FDA', 'NCI_GAIA', 'HL7V3.0', 'NEU', 'PSY', 'SPN', 'AIR', 'GO', 'CCC', 'UMD', 
    'NIC', 'ALT', 'NCI_EDQM-HC', 'NCI_INC', 'LNC', 'JABL', 'NUCCPT', 'ICF', 'ICF-CY', 'NCI_BRIDG_3_0_3', 
    'NCI_BRIDG_5_3', 'NANDA-I', 'PNDS', 'NOC', 'OMS', 'NCI_CTEP-SDC', 'NCI_DICOM', 'NCI_KEGG', 
    'NCI_BioC', 'MCM', 'SOP', 'NCI_CTCAE_5', 'NCI_CTCAE_3', 'MDR'
]

# Define columns for MRDEF.RRF based on MRDEF.ctl
COLUMNS_DEF = ['CUI', 'AUI', 'ATUI', 'SATUI', 'SAB', 'DEF', 'SUPPRESS', 'CVF']

# Define columns for MRCONSO.RRF based on the new data structure
COLUMNS_CONSO = [
    'cui', 'language', 'term_status', 'lui', 'string_type', 'string_identifier', 'is_preferred',
    'aui', 'source_aui', 'source_cui', 'source_descriptor_dui', 'source', 
    'source_term_type', 'source_code', 'source_name', 'srl', 'suppress', 'cvf'
]

def extract_cuis_from_mrsty(input_file):
    disorder_cui = set()
    disease_cui = set()
    
    with open(input_file, 'r') as infile:
        for line in infile:
            fields = line.strip().split('|')
            if len(fields) >= len(COLUMNS_STY):
                cui, sty = fields[0], fields[3]
                if sty in DISORDER_ENTITIES:
                    disorder_cui.add(cui)
                if sty == 'Disease or Syndrome':
                    disease_cui.add(cui)
    return disorder_cui, disease_cui

def process_mrdef(input_file, output_file_disorder, output_file_disease, disorder_cui, disease_cui):
    with open(input_file, 'r') as infile, \
         open(output_file_disorder, 'w', newline='', encoding='utf-8') as outfile_disorder, \
         open(output_file_disease, 'w', newline='', encoding='utf-8') as outfile_disease:
        
        writer_disorder = csv.writer(outfile_disorder, delimiter='\t')
        writer_disease = csv.writer(outfile_disease, delimiter='\t')
        
        # Write headers
        writer_disorder.writerow(['CUI', 'SAB', 'DEF'])
        writer_disease.writerow(['CUI', 'SAB', 'DEF'])
        
        for line in infile:
            fields = line.strip().split('|')
            if len(fields) >= len(COLUMNS_DEF):
                cui, sab, definition = fields[0], fields[4], fields[5]
                if sab in ENGLISH_SOURCES:
                    if cui in disorder_cui:
                        writer_disorder.writerow([cui, sab, definition])
                    if cui in disease_cui:
                        writer_disease.writerow([cui, sab, definition])

def process_mrconso(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(COLUMNS_CONSO)  # Write the header to the output file
        
        for line in infile:
            fields = line.strip().split('|')
            if len(fields) >= len(COLUMNS_CONSO) and fields[1].strip() == 'ENG':  # Check if the line is in English and has at least the expected number of columns
                writer.writerow(fields[:len(COLUMNS_CONSO)])  # Write only the expected number of columns

def main():
    disorder_cui, disease_cui = extract_cuis_from_mrsty(INPUT_FILE_STY)
    process_mrdef(INPUT_FILE_DEF, OUTPUT_FILE_DISORDER, OUTPUT_FILE_DISEASE, disorder_cui, disease_cui)
    process_mrconso(INPUT_FILE_CONSO, OUTPUT_FILE_CONSO)
    print(f"Processing complete. Saved to {OUTPUT_FILE_DISORDER}, {OUTPUT_FILE_DISEASE}, and {OUTPUT_FILE_CONSO}")

if __name__ == '__main__':
    main()
