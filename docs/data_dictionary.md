## Data Dictionary

### Table: `patient`

| Field Name | Data Type | Key / Constraint | Description |
|---|---|---|---|
| `patient_id` | `VARCHAR(100)` | Primary Key, Not Null | Unique identifier for each patient from `PATIENT_ID`. |
| `gender` | `VARCHAR(10)` | Nullable if unknown | Patient's recorded gender. |
| `race_id` | `INT` | Foreign Key → `race(race_id)` | Links the patient to a race category. |
| `ethnicity_id` | `INT` | Foreign Key → `ethnicity(ethnicity_id)` | Links the patient to an ethnicity category. |

### Table: `race`

| Field Name | Data Type | Key / Constraint | Description |
|---|---|---|---|
| `race_id` | `INT` | Primary Key, Not Null | Unique identifier for each race category. |
| `race` | `VARCHAR(100)` | Unique, Not Null | Race category from the data file. |

### Table: `ethnicity`

| Field Name | Data Type | Key / Constraint | Description |
|---|---|---|---|
| `ethnicity_id` | `INT` | Primary Key, Not Null | Unique identifier for each ethnicity category. |
| `ethnicity` | `VARCHAR(100)` | Unique, Not Null | Ethnicity category from the data file. |

### Table: `sample`

| Field Name | Data Type | Key / Constraint | Description |
|---|---|---|---|
| `sample_id` | `VARCHAR(100)` | Primary Key, Not Null | Unique identifier for each tumor sample. |
| `patient_id` | `VARCHAR(100)` | Foreign Key → `patient(patient_id)`, Not Null | Links the sample to the patient. |
| `disease_id` | `INT` | Foreign Key → `disease(disease_id)`, Not Null | Links the sample to a cancer type. |
| `age` | `INT` | Nullable if unknown | Patient age at the time of sample collection or record. |
| `stage_id` | `INT` | Foreign Key → `stage(stage_id)` | Links the sample to a cancer stage. |
| `tumor_tissue` | `VARCHAR(150)` | Nullable if unknown | Tumor source description. |

### Table: `disease`

| Field Name | Data Type | Key / Constraint | Description |
|---|---|---|---|
| `disease_id` | `INT` | Primary Key, Not Null | Unique identifier for each disease or cancer type. |
| `cancer_type` | `VARCHAR(100)` | Unique, Not Null | Standardized cancer type name. |

### Table: `stage`

| Field Name | Data Type | Key / Constraint | Description |
|---|---|---|---|
| `stage_id` | `INT` | Primary Key, Not Null | Unique identifier for each cancer stage. |
| `stage` | `VARCHAR(20)` | Unique, Not Null | Cancer stage label such as I, II, III, IIIA, IIIB, IIIC, IV, IVA, or IVB. |

### Table: `gene`

| Field Name | Data Type | Key / Constraint | Description |
|---|---|---|---|
| `gene_id` | `INT` | Primary Key, Not Null | Internal unique identifier for each gene. |
| `hugo_symbol` | `VARCHAR(100)` | Unique but nullable if unknown | Standard HUGO gene symbol, such as `TP53`, `BRCA1`, or `EGFR`. |
| `entrez_id` | `INT` | Unique but nullable if unknown | External Entrez gene identifier. |
| `chromosome` | `VARCHAR(10)` | Nullable | Chromosome where the gene is located. |

### Table: `mutation`

| Field Name | Data Type | Key / Constraint | Description |
|---|---|---|---|
| `mutation_id` | `INT` | Primary Key, Not Null | Unique identifier for each mutation record. |
| `gene_id` | `INT` | Foreign Key → `gene(gene_id)`, Not Null | Links the mutation to the gene where it occurs. |
| `HGVSp` | `VARCHAR(100)` | Nullable | Protein-level mutation notation. |
| `variantClassification` | `VARCHAR(100)` | Nullable | Functional classification of the mutation. |
| `codon_change` | `VARCHAR(100)` | Nullable | Codon-level change caused by the mutation. |
| `start_pos` | `INT` | Nullable | Genomic start position of the mutation. |
| `end_pos` | `INT` | Nullable | Genomic end position of the mutation. |
| `reference_allele` | `TEXT` | Nullable | Reference genome allele. |
| `tumor_allele1` | `TEXT` | Nullable | First observed tumor allele. |
| `tumor_allele2` | `TEXT` | Nullable | Second observed tumor allele. |
| `polyphen_score` | `VARCHAR(100)` | Nullable | Predicted effect score of the mutation. |
| `varian_type` | `VARCHAR(10)` | Nullable | General type of mutation or variant. |
| `dbSNP` | `VARCHAR(100)` | Nullable | dbSNP reference identifier. |

### Table: `sample_mutation`

| Field Name | Data Type | Key / Constraint | Description |
|---|---|---|---|
| `sample_mutation_id` | `INT` | Primary Key, Not Null | Unique identifier for each sample-mutation relationship. |
| `sample_id` | `VARCHAR(100)` | Foreign Key → `sample(sample_id)`, Not Null | Links the record to a sample. |
| `mutation_id` | `INT` | Foreign Key → `mutation(mutation_id)`, Not Null | Links the record to a mutation. |

### Table: `expression`

| Field Name | Data Type | Key / Constraint | Description |
|---|---|---|---|
| `sample_id` | `VARCHAR(100)` | Composite Primary Key, Foreign Key → `sample(sample_id)` | Identifies the sample for the expression measurement. |
| `gene_id` | `INT` | Composite Primary Key, Foreign Key → `gene(gene_id)` | Identifies the gene for the expression measurement. |
| `expressionValue` | `DECIMAL(20,8)` | Not Null | Measured expression value for a specific gene in a specific sample. |
