# Project Title

## Project Summary
This project builds a relational database using TCGA hepatocellular carcinoma clinical and genomic data downloaded from cBioPortal. The database is designed so users can search by patient ID, sample ID, gene, mutation, or SNP identifier to retrieve clinical information, mutation details, and raw mRNA expression values

## Tools and Technologies

- **DBMS:** MySQL / MariaDB
- **Database Interface:** phpMyAdmin
- **Programming Language:** Python 3.13
- **Python Libraries:** pandas
- **IDE:** PyCharm 2025.2.2
- **Virtualization / Local Server:** Oracle VirtualBox 7.2.0 or XAMPP
- **Data Source:** cBioPortal / TCGA Liver Hepatocellular Carcinoma dataset
- **Version Control:** Git and GitHub

## Repository Structure
```text
.
├── data/
│   ├── raw/                    # Original TCGA/cBioPortal data files. Inside is the link to Figshare due to large data size
│   └── cleaned/                # Cleaned and mapped TSV files. Inside is the link to Figshare due to large data size
├── diagrams/
│   ├── ER_diagram.png               # Entity-relationship diagram of the database
│   ├── README.md                 
│   ├── database_design.png          # Complete database design diagram
│   └── loaded_database.png          # Screenshot of populated database tables
├── scripts/
│   ├── 01_clean_data.py             # Cleans raw files and standardizes missing values
│   ├── 02_map_genes.py              # Maps outdated or missing HUGO symbols and Entrez IDs
│   ├── 03_fix_mutation_from_expression.py
│   │                                # Fixes remaining mutation Entrez IDs using expression data
│   └── 04_write_sql.py              # Generates SQL INSERT statements from cleaned data
├── sql/                            # All sql files are uploaded to figshare due to large data size. Link is found in README
│   └── README.md    #contain link to figshare to download 01_create_schema.sql, 02_load_cleaned_data.sql, database_dump.sql
├── docs/
│   └── Project_Documentation.pdf    # Final project documentation
├── README.md                        # Project overview and data dictionary
└── .gitignore                       # Files and folders ignored by Git
```

## Data Sources
Data was from TCGA hepatocellular carcinoma clinical and genomic data downloaded from cBioPortal.

## How to Recreate the Database
1. Install Oracle VirtualBox, local XAMPP, MySQL, Python 3.13, and an IDE such as PyCharm.
2. Download the raw data files from the GitHub repository or from cBioPortal.
3. Run `01_clean_data.py` separately on:
   - `data_clinical_patient.txt`
   - `data_clinical_sample.txt`
   - `data_mutations.txt`
   - `data_mrna_seq_v2_rsem.txt`
4. Run `02_map_genes.py` using `hgnc_complete_set.txt` as the reference file.
5. Run `03_fix_mutation_from_expression.py` to finalize mutation gene identifier mapping.
6. Run `04_write_sql.py` to generate `02_load_cleaned_data.sql`.
7. Run `01_create_schema.sql` to create the database schema.
8. Run `02_load_cleaned_data.sql` to populate the database.
9. Access the database through phpMyAdmin or MySQL.

Example MySQL command:

```bash
mysql -h 127.0.0.1 -P 3306 -u root -pinstructor
```
## Documentation and Diagrams
The write-up, data dictionary, script execution order, decisions and limitation can be found in the directory docs. Diagrams can be found in directory diagrams.
