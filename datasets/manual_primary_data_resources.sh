# DRUGBANK
# Required DrugBank CSV files stored in PostgreSQL database called db
# Database: DrugBank, Output: 12 drug feature files, drug_protein.csv
echo "Parsing DrugBank data..."
python processing_scripts/drugbank_parsedb.py
python processing_scripts/drugbank_drug_protein.py


# DRUG CENTRAL
# Database: Drug Central, Output: drug_disease.tsv, structures.tsv, structure_type.tsv, dc_features.tsv
# https://unmtid-shinyapps.net/download/drugcentral.dump.11012023.sql.gz
curl "https://unmtid-shinyapps.net/download/drugcentral.dump.11012023.sql.gz" -o data/input/drugcentral/drugcentral.dump.11012023.sql.gz
gunzip data/input/drugcentral/drugcentral.dump.11012023.sql.gz
dropdb drugcentral
createdb drugcentral
psql drugcentral < data/input/drugcentral/drugcentral.dump.11012023.sql
psql -d drugcentral -c "SELECT DISTINCT * FROM structures RIGHT JOIN (SELECT * FROM omop_relationship WHERE relationship_name IN ('indication', 'contraindication', 'off-label use')) AS drug_disease ON structures.id = drug_disease.struct_id;" -P format=unaligned -P fieldsep=$'\t' -o data/output/drugcentral/drug_disease.tsv
psql -d drugcentral -c "COPY (SELECT * FROM structures) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > data/input/drugcentral/structures.tsv
psql -d drugcentral -c "COPY (SELECT * FROM structure_type) TO STDOUT WITH CSV HEADER DELIMITER E'\t'" > data/input/drugcentral/structure_type.tsv
python processing_scripts/drugcentral_feature.py


# Database: UMLS, Script: umls.py, map_umls_mondo.py, Output: umls.tsv, umls_def_disorder.tsv, umls_def_disease.tsv, umls_mondo.tsv
# The UMLS Metathesaurus will need to be downloaded manually after authentication at uts.nlm.nih.gov.
mv ~/Downloads/umls-2024AA-metathesaurus-full.zip data/input/umls/
unzip data/input/umls/umls-2024AA-metathesaurus-full.zip -d data/input/umls/
python processing_scripts/umls.py
python processing_scripts/map_umls_mondo.py