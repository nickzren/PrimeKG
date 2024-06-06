import os
import pandas as pd
from sqlalchemy import create_engine

# Constants
OUTPUT_DIR = 'data/output/drugbank/'

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Function to fetch data from PostgreSQL and save to TSV
def fetch_and_save_to_tsv(query, output_file, engine):
    try:
        # Execute the query and fetch data
        df = pd.read_sql_query(query, engine)
        
        # Drop any rows with missing values and remove duplicates
        df.dropna(axis=0, how='any', inplace=True)
        df.drop_duplicates(inplace=True)
        
        # Save the DataFrame to a TSV file
        df.to_csv(output_file, sep='\t', index=False)
        
        print(f"Data saved to {output_file}")
    
    except Exception as e:
        print(f"Error: {e}")

# Queries to fetch data
queries = {
    "volume_of_distribution": "SELECT drugbank_id AS ID, volume_of_distribution FROM db.drugs JOIN db.pharmacologies ON db.drugs.id = db.pharmacologies.drug_id WHERE volume_of_distribution IS NOT NULL",
    "half_life": "SELECT drugbank_id AS ID, half_life FROM db.drugs JOIN db.pharmacologies ON db.drugs.id = db.pharmacologies.drug_id WHERE half_life IS NOT NULL",
    "pharmacodynamics": "SELECT drugbank_id AS ID, pharmacodynamics FROM db.drugs JOIN db.pharmacologies ON db.drugs.id = db.pharmacologies.drug_id WHERE pharmacodynamics IS NOT NULL",
    "groups": """
        SELECT drugbank_id AS ID, 'experimental' AS group FROM db.drugs WHERE experimental IS TRUE
        UNION ALL
        SELECT drugbank_id AS ID, 'investigational' AS group FROM db.drugs WHERE investigational IS TRUE
        UNION ALL
        SELECT drugbank_id AS ID, 'approved' AS group FROM db.drugs WHERE approved IS TRUE
        UNION ALL
        SELECT drugbank_id AS ID, 'vet_approved' AS group FROM db.drugs WHERE vet_approved IS TRUE
        UNION ALL
        SELECT drugbank_id AS ID, 'withdrawn' AS group FROM db.drugs WHERE withdrawn IS TRUE
        UNION ALL
        SELECT drugbank_id AS ID, 'illicit' AS group FROM db.drugs WHERE illicit IS TRUE
        UNION ALL
        SELECT drugbank_id AS ID, 'nutraceutical' AS group FROM db.drugs WHERE nutraceutical IS TRUE
    """,
    "description": "SELECT drugbank_id AS ID, description FROM db.drugs WHERE description IS NOT NULL",
    "state": "SELECT drugbank_id AS ID, state FROM db.drugs WHERE state IS NOT NULL",
    "indication": "SELECT drugbank_id AS ID, indication FROM db.drugs JOIN db.pharmacologies ON db.drugs.id = db.pharmacologies.drug_id WHERE indication IS NOT NULL",
    "protein_binding": "SELECT drugbank_id AS ID, protein_binding FROM db.drugs JOIN db.pharmacologies ON db.drugs.id = db.pharmacologies.drug_id WHERE protein_binding IS NOT NULL",
    "pathways": """
        SELECT d.drugbank_id AS ID, p.name || '|' || p.category AS pathway
        FROM db.drugs d
        JOIN db.pathway_drugs pd ON d.id = pd.drug_id
        JOIN db.pathways p ON pd.smpdb_id = p.smpdb_id
    """,
    "mechanism_of_action": "SELECT drugbank_id AS ID, mechanism_of_action FROM db.drugs JOIN db.pharmacologies ON db.drugs.id = db.pharmacologies.drug_id WHERE mechanism_of_action IS NOT NULL",
    "chemical_taxonomy": """
        SELECT d.drugbank_id AS ID, 'Kingdom' AS Class, t.kingdom AS Chem FROM db.drugs d
        JOIN db.taxonomy t ON d.id = t.drug_id
        UNION ALL
        SELECT d.drugbank_id AS ID, 'Superclass' AS Class, t.superklass AS Chem FROM db.drugs d
        JOIN db.taxonomy t ON d.id = t.drug_id
        UNION ALL
        SELECT d.drugbank_id AS ID, 'Class' AS Class, t.klass AS Chem FROM db.drugs d
        JOIN db.taxonomy t ON d.id = t.drug_id
        UNION ALL
        SELECT d.drugbank_id AS ID, 'Subclass' AS Class, t.subklass AS Chem FROM db.drugs d
        JOIN db.taxonomy t ON d.id = t.drug_id
        UNION ALL
        SELECT d.drugbank_id AS ID, 'Direct parent' AS Class, t.direct_parent AS Chem FROM db.drugs d
        JOIN db.taxonomy t ON d.id = t.drug_id
        ORDER BY ID, Class
    """,
    "drug_categories": """
        SELECT d.drugbank_id AS ID, c.title AS Drugcat
        FROM db.drug_categorizations dc
        JOIN db.categories c ON dc.category_id = c.id
        JOIN db.drugs d ON dc.drug_id = d.id
    """,
    "atc_codes": """
        SELECT DISTINCT d.drugbank_id AS ID, 
            STRING_AGG(DISTINCT CASE WHEN cm.vocabulary_level = 1 THEN cm.title ELSE NULL END, '|') AS Level1,
            STRING_AGG(DISTINCT CASE WHEN cm.vocabulary_level = 2 THEN cm.title ELSE NULL END, '|') AS Level2,
            STRING_AGG(DISTINCT CASE WHEN cm.vocabulary_level = 3 THEN cm.title ELSE NULL END, '|') AS Level3,
            STRING_AGG(DISTINCT CASE WHEN cm.vocabulary_level = 4 THEN cm.title ELSE NULL END, '|') AS Level4,
            STRING_AGG(DISTINCT CASE WHEN cm.vocabulary_level = 5 THEN cm.title ELSE NULL END, '|') AS Level5
        FROM db.drugs d
        LEFT JOIN db.drug_categorizations dc ON dc.drug_id = d.id
        LEFT JOIN db.category_mappings cm ON cm.category_id = dc.category_id
        WHERE cm.vocabulary = 'ATC'
        GROUP BY d.drugbank_id 
        ORDER BY d.drugbank_id
    """,
    "drugbank_vocabulary": """
        SELECT
            d.drugbank_id,
            STRING_AGG(DISTINCT an.number, ',') AS AccessionNumbers,
            d.name, 
            d.cas_number, 
            STRING_AGG(DISTINCT u.unii, ',') AS UNII,
            STRING_AGG(DISTINCT ds.synonym, ',') AS Synonyms,
            d.moldb_inchikey AS StandardInChiKey
        FROM db.drugs d
        JOIN db.accession_numbers an ON an.record_id = d.id
        JOIN db.drug_synonyms ds ON ds.drug_id = d.id
        JOIN db.uniis u ON u.record_id = d.id
        WHERE an.record_type = 'Drug' AND u.record_type = 'Drug'
        GROUP BY d.drugbank_id, d.name, d.cas_number, d.moldb_inchikey
        ORDER BY d.drugbank_id
    """,
    "drug_bonds": """
        SELECT
            b.id,
            p.name, 
            p.uniprot_id, 
            be.organism AS species, 
            STRING_AGG(DISTINCT d.drugbank_id, '; ') AS "Drug IDs",
            b.type AS bond_type
        FROM db.drugs d
        JOIN db.bonds b ON b.drug_id = d.id
        JOIN db.bio_entities be ON be.biodb_id = b.biodb_id
        JOIN db.bio_entity_components bec ON bec.biodb_id = b.biodb_id
        JOIN db.polypeptides p ON p.uniprot_id = bec.component_id
        WHERE be.organism = 'Humans' AND b.type IN ('CarrierBond', 'EnzymeBond', 'TargetBond', 'TransporterBond')
        GROUP BY b.id, p.name, p.uniprot_id, be.organism, b.type
        ORDER BY b.id
    """,
    "drugs": """
        SELECT drugbank_id AS ID, type AS drug_type, cas_number FROM db.drugs
    """
}

def main():
    try:
        # Retrieve connection string from environment variable
        connection_string = os.getenv('DB_CONNECTION_STRING')
        if not connection_string:
            raise ValueError("DB_CONNECTION_STRING environment variable not set")

        # Ensure using the correct dialect
        connection_string = connection_string.replace('postgres://', 'postgresql+psycopg2://')
        
        # Create SQLAlchemy engine
        engine = create_engine(connection_string)
        
        # Execute the queries and save the results to TSV files
        for key, query in queries.items():
            output_file = os.path.join(OUTPUT_DIR, f"{key}.tsv")
            fetch_and_save_to_tsv(query, output_file, engine)

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
