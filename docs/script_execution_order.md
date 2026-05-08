## Scripts execution order

### Python Scripts

| Order | Script | Description | Output |
|---|---|---|---|
| 1 | `01_clean_data.py` | Performs initial cleaning, checks duplicates, replaces missing-value placeholders, converts raw `.txt` files to `.tsv`, and removes empty columns. | Cleaned TSV files |
| 2 | `02_map_genes.py` | Updates outdated or missing HUGO symbols and Entrez IDs using `hgnc_complete_set.txt`. | `expression_mapped_final.tsv`, `mutation_mapped.tsv` |
| 3 | `03_fix_mutation_from_expression.py` | Uses the cleaned expression file to fix remaining invalid mutation Entrez IDs where possible. | `mutation_mapped_fixed_final.tsv` |
| 4 | `04_write_sql.py` | Maps cleaned data into database schema fields and generates SQL insert statements. | `02_load_cleaned_data.sql` |

### SQL Scripts

| Order | Script | Description |
|---|---|---|
| 1 | `01_create_schema.sql` | Creates the database and tables. |
| 2 | `02_load_cleaned_data.sql` | Loads cleaned data into the database. |
| N/A | `database_dump` | Backup or transfer script for recreating the database. |
