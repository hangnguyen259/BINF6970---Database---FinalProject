## Source Column Mapping

### Patient Data

| Source Column | Database Field |
|---|---|
| `PATIENT_ID` | `patient_id` |
| `SEX` | `gender` |
| `RACE` | `race` |
| `ETHNICITY` | `ethnicity` |
| `AGE` | `age` |
| `CANCER_TYPE_ACRONYM` | `cancer_type` |
| `AJCC_PATHOLOGIC_TUMOR_STAGE` | `stage` |

### Sample Data

| Source Column | Database Field |
|---|---|
| `PATIENT_ID` | `patient_id` |
| `SAMPLE_ID` | `sample_id` |
| `TUMOR_TISSUE_SITE` | `tumor_tissue` |

### Mutation Data

| Source Column | Database Field |
|---|---|
| `HUGO_SYMBOL` | `hugo_symbol` |
| `ENTREZ_GENE_ID` | `entrez_id` |
| `CHROMOSOME` | `chromosome` |
| `START_POSITION` | `start_pos` |
| `END_POSITION` | `end_pos` |
| `VARIANT_CLASSIFICATION` | `variantClassification` |
| `VARIANT_TYPE` | `varian_type` |
| `REFERENCE_ALLELE` | `reference_allele` |
| `TUMOR_SEQ_ALLELE1` | `tumor_allele1` |
| `TUMOR_SEQ_ALLELE2` | `tumor_allele2` |
| `dbSNP_RS` | `dbSNP` |
| `TUMOR_SAMPLE_BARCODE` | `sample_id` |
| `HGVSP_SHORT` | `HGVSp` |
| `CODONS` | `codon_change` |
| `POLYPHEN` | `polyphen_score` |

### Expression Data

| Source Column | Database Field |
|---|---|
| `HUGO_symbol` | `hugo_symbol` |
| `Entrez_ID` | `entrez_id` |
| Tumor sample barcode columns | Expression values by gene and sample |

## Data Cleaning Summary

The data cleaning workflow standardizes missing values, removes empty columns, checks for duplicate records, and updates outdated or missing gene identifiers.

Key cleaning decisions:

- Missing-value placeholders such as `.` were replaced with standardized null values.
- Mutation and expression gene identifiers were mapped using `hgnc_complete_set.txt`.
- Invalid Entrez IDs in mutation data were corrected where possible using the expression file as a reference.
- Mutation rows were not dropped if they had at least one valid gene identifier.
- Raw expression values, including zero values, were kept to preserve flexibility for downstream analysis.
  ## Limitations

Current limitations include:

- The database supports one disease per sample, which works for the liver cancer dataset but may need redesign for multi-disease datasets.
- `sample_mutation` uses a surrogate primary key instead of a composite key.
- Some Python scripts contain hardcoded file names.
- Expression data are raw values and may require normalization before downstream analysis.
- Only a subset of available TCGA/cBioPortal files was included in the database.
