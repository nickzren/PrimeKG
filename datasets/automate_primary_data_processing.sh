#!/bin/bash

LOG_FILE="automate_primary_data_processing.log"

# Recreate the log file every new run of the script
echo "Starting new log file at $(date)" > "$LOG_FILE"

# Redirect stdout and stderr to the log file
exec > >(tee -a "$LOG_FILE") 2>&1

# Create the main directories, including subdirectories for input and output
echo "Making required directories..."
dirs=("bgee" "ctd" "disgenet" "drugcentral" "drugbank" "vocab" "ncbigene" "go" "hpo" "mondo" "reactome" "sider" "uberon" "umls")
for dir in "${dirs[@]}"; do
  mkdir -p "data/input/$dir" "data/output/$dir"
done

# GENE NAMES
echo "Downloading gene names..."
curl "https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=gd_app_name&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_pub_eg_id&col=md_eg_id&col=md_prot_id&col=md_mim_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit" -o data/input/vocab/gene_names.csv
curl "https://www.genenames.org/cgi-bin/download/custom?col=md_eg_id&col=md_prot_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit" -o data/input/vocab/gene_map.csv

# BGEE
echo "Downloading and processing BGEE data..."
curl https://www.bgee.org/ftp/current/download/calls/expr_calls/Homo_sapiens_expr_advanced.tsv.gz -o data/input/bgee/Homo_sapiens_expr_advanced.tsv.gz
gunzip data/input/bgee/Homo_sapiens_expr_advanced.tsv.gz
python processing_scripts/bgee.py

# COMPARATIVE TOXICOGENOMICS DATABASE
echo "Downloading and processing Comparative Toxicogenomics Database data..."
curl https://ctdbase.org/reports/CTD_exposure_events.csv.gz -o data/input/ctd/CTD_exposure_events.csv.gz
gunzip data/input/ctd/CTD_exposure_events.csv.gz
python processing_scripts/ctd.py

# DISGENET
echo "Downloading and processing DisGeNET data..."
curl https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz -o data/input/disgenet/curated_gene_disease_associations.tsv.gz
gunzip data/input/disgenet/curated_gene_disease_associations.tsv.gz

# ENTREZ GENE
echo "Downloading Entrez Gene data..."
curl https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz -o data/input/ncbigene/gene2go.gz
gunzip data/input/ncbigene/gene2go.gz
python processing_scripts/ncbigene.py

# Gene Ontology
echo "Downloading and processing Gene Ontology data..."
curl -L http://purl.obolibrary.org/obo/go/go-basic.obo -o data/input/go/go-basic.obo
python processing_scripts/go.py

# Human Phenotype Ontology
echo "Downloading and processing Human Phenotype Ontology data..."
curl -L http://purl.obolibrary.org/obo/hp.obo -o data/input/hpo/hp.obo
python processing_scripts/hpo.py

curl -L http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa -o data/input/hpo/phenotype.hpoa
python processing_scripts/hpoa.py

# MONDO
echo "Downloading and processing MONDO data..."
curl -L http://purl.obolibrary.org/obo/MONDO.obo -o data/input/mondo/mondo.obo
python processing_scripts/mondo.py

# Reactome
echo "Downloading and processing Reactome data..."
curl https://reactome.org/download/current/ReactomePathways.txt -o data/input/reactome/ReactomePathways.txt
curl https://reactome.org/download/current/ReactomePathwaysRelation.txt -o data/input/reactome/ReactomePathwaysRelation.txt
curl https://reactome.org/download/current/NCBI2Reactome.txt -o data/input/reactome/NCBI2Reactome.txt
python processing_scripts/reactome.py

# SIDER
echo "Downloading and processing SIDER data..."
curl http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz -o data/input/sider/meddra_all_se.tsv.gz
gunzip data/input/sider/meddra_all_se.tsv.gz
curl http://sideeffects.embl.de/media/download/drug_atc.tsv -o data/input/sider/drug_atc.tsv
python processing_scripts/sider.py

# UBERON
echo "Downloading and processing UBERON data..."
curl -L http://purl.obolibrary.org/obo/uberon/ext.obo -o data/input/uberon/ext.obo
python processing_scripts/uberon.py

# Create directories for output files and PPI data
echo "Setting up directories for knowledge graph and PPI data..."
mkdir -p data/output/kg/auxillary data/output/ppi
