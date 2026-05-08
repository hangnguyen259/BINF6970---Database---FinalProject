import re
from pathlib import Path
import pandas as pd

# =========================
# CONFIG
# =========================

PATIENT_FILE = "patient_data_cleaned_final.tsv"
SAMPLE_FILE = "sample_data_cleaned_final.tsv"
MUTATION_FILE = "mutation_mapped_fixed_final.tsv"
EXPRESSION_FILE = "expression_mapped_final.tsv"

OUTPUT_SQL = "load_w_ids2.sql"
DB_NAME = "FinalProject"

SEP = "\t"
BATCH_SIZE = 500


# =========================
# HELPERS
# =========================

def clean_colname(col):
    """Clean header formatting without changing capitalization."""
    col = str(col).strip()
    col = col.replace(" ", "_")
    col = re.sub(r"[^\w]+", "_", col)
    col = re.sub(r"_+", "_", col)
    return col.strip("_")


def header_key(col):
    """Case-insensitive key used only for recognizing known input columns."""
    return clean_colname(col).lower()


def assert_no_duplicate_columns(df, name):
    """
    Fail early if cleaned headers create duplicate column names.

    This avoids pandas returning a DataFrame from df[c] when duplicate
    column names exist, which can cause pd.isna(...) ambiguity errors.
    """
    dupes = df.columns[df.columns.duplicated()].tolist()
    if dupes:
        raise ValueError(
            f"{name} file has duplicate columns after header cleaning: {dupes}. "
            "Rename or remove the duplicate columns before running this script."
        )


def norm_text(x):
    """Normalize a single cell value to None or stripped string."""
    if x is None:
        return None
    if pd.isna(x):
        return None
    s = str(x).strip()
    return s if s != "" else None


def norm_key(x):
    """Normalize IDs/keys to None or stripped string."""
    if x is None:
        return None
    if pd.isna(x):
        return None
    s = str(x).strip()
    return s if s != "" else None


def norm_gene_symbol(x):
    """Normalize gene symbols so TP53, tp53, and ' TP53 ' match."""
    s = norm_key(x)
    return s.upper() if s else None


def norm_entrez_id(x):
    """
    Normalize Entrez IDs so values like 7157, '7157', and '7157.0' match.
    """
    s = norm_key(x)
    if not s:
        return None

    try:
        f = float(s)
        if f.is_integer():
            return str(int(f))
    except ValueError:
        pass

    return s


def sql_quote(val):
    if val is None or pd.isna(val):
        return "NULL"
    s = str(val).replace("\\", "\\\\").replace("'", "''")
    return f"'{s}'"


def sql_number(val):
    if val is None or pd.isna(val) or str(val).strip() == "":
        return "NULL"
    return str(val)


def make_id_map(values, start_id=1):
    return {value: i for i, value in enumerate(values, start=start_id)}


def batched(seq, size):
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def multirow_insert(table, columns, rows):
    if not rows:
        return []

    statements = []
    col_sql = ", ".join(columns)

    for batch in batched(rows, BATCH_SIZE):
        values_sql = ",\n".join(
            "(" + ", ".join(row) + ")" for row in batch
        )
        statements.append(
            f"INSERT INTO {table} ({col_sql}) VALUES\n{values_sql};"
        )

    return statements


def sort_key_tuple(t):
    return tuple("" if v is None or pd.isna(v) else str(v) for v in t)


def is_tcga_sample_column(col):
    """
    Return True if an expression file column looks like a TCGA sample barcode.

    Header cleaning converts TCGA-XX-XXXX... to TCGA_XX_XXXX..., so this
    accepts either hyphen or underscore separators.
    """
    s = str(col).strip()
    return re.match(r"^TCGA[-_][A-Z0-9]{2}[-_][A-Z0-9]{4}", s, re.IGNORECASE) is not None


# =========================
# READ FILES
# =========================

patient_df = pd.read_csv(PATIENT_FILE, sep=SEP, dtype=str)
sample_df = pd.read_csv(SAMPLE_FILE, sep=SEP, dtype=str)
mutation_df = pd.read_csv(MUTATION_FILE, sep=SEP, dtype=str)
expression_df = pd.read_csv(EXPRESSION_FILE, sep=SEP, dtype=str)

patient_df.columns = [clean_colname(c) for c in patient_df.columns]
sample_df.columns = [clean_colname(c) for c in sample_df.columns]
mutation_df.columns = [clean_colname(c) for c in mutation_df.columns]
expression_df.columns = [clean_colname(c) for c in expression_df.columns]

assert_no_duplicate_columns(patient_df, "patient")
assert_no_duplicate_columns(sample_df, "sample")
assert_no_duplicate_columns(mutation_df, "mutation")
assert_no_duplicate_columns(expression_df, "expression")

print("Patient headers:", list(patient_df.columns))
print("Sample headers:", list(sample_df.columns))
print("Mutation headers:", list(mutation_df.columns))
print("Expression headers first 10:", list(expression_df.columns[:10]))


# =========================
# HEADER MAPPING
# =========================

patient_map = {}
for c in patient_df.columns:
    k = header_key(c)
    if k == "patient_id":
        patient_map[c] = "patient_id"
    elif k in ("sex"):
        patient_map[c] = "gender"
    elif k == "race":
        patient_map[c] = "race"
    elif k == "ethnicity":
        patient_map[c] = "ethnicity"
    elif k == "age":
        patient_map[c] = "age"
    elif k in ("cancer_type_acronym"):
        patient_map[c] = "cancer_type"
    elif k in ("ajcc_pathologic_tumor_stage"):
        patient_map[c] = "stage"

patient_df = patient_df.rename(columns=patient_map)
#patient_df = coalesce_duplicate_columns(patient_df, "patient")

for col in ["patient_id", "gender", "race", "ethnicity"]:
    if col not in patient_df.columns:
        raise ValueError(f"Missing expected patient column: {col}")

for optional_col in ["stage", "age", "cancer_type"]:
    if optional_col not in patient_df.columns:
        patient_df[optional_col] = None


sample_map = {}
for c in sample_df.columns:
    k = header_key(c)
    if k == "patient_id":
        sample_map[c] = "patient_id"
    elif k == "sample_id":
        sample_map[c] = "sample_id"
    elif k in ("tumor_tissue_site"):
        sample_map[c] = "tumor_tissue"

sample_df = sample_df.rename(columns=sample_map)
#sample_df = coalesce_duplicate_columns(sample_df, "sample")

for col in ["patient_id", "sample_id", "tumor_tissue"]:
    if col not in sample_df.columns:
        raise ValueError(f"Missing expected sample column: {col}")


mutation_map = {}
for c in mutation_df.columns:
    k = header_key(c)
    if k == "hugo_symbol":
        mutation_map[c] = "hugo_symbol"
    elif k == "entrez_gene_id":
        mutation_map[c] = "entrez_id"
    elif k == "chromosome":
        mutation_map[c] = "chromosome"
    elif k == "start_position":
        mutation_map[c] = "start_pos"
    elif k == "end_position":
        mutation_map[c] = "end_pos"
    elif k == "variant_classification":
        mutation_map[c] = "variant_classification"
    elif k == "variant_type":
        mutation_map[c] = "variant_type"
    elif k == "reference_allele":
        mutation_map[c] = "reference_allele"
    elif k == "tumor_seq_allele1":
        mutation_map[c] = "tumor_allele1"
    elif k == "tumor_seq_allele2":
        mutation_map[c] = "tumor_allele2"
    elif k in ("dbsnp_rs"):
        mutation_map[c] = "dbsnp"
    elif k in ("tumor_sample_barcode"):
        mutation_map[c] = "sample_id"
    elif k == "hgvsp_short":
        # Use only HGVSp_Short for the SQL HGVSp field.
        # Do not map a general HGVSp column, because it may contain longer protein annotations.
        mutation_map[c] = "hgvsp"
    elif k in ("codons"):
        mutation_map[c] = "codon_change"
    elif k in ("polyphen"):
        mutation_map[c] = "polyphen_score"

mutation_df = mutation_df.rename(columns=mutation_map)
#mutation_df = coalesce_duplicate_columns(mutation_df, "mutation")

required_mut_cols = [
    "hugo_symbol", "entrez_id", "chromosome", "start_pos", "end_pos",
    "variant_classification", "variant_type", "reference_allele",
    "tumor_allele1", "tumor_allele2", "dbsnp", "sample_id",
    "hgvsp", "codon_change", "polyphen_score"
]

for col in required_mut_cols:
    if col not in mutation_df.columns:
        raise ValueError(f"Missing expected mutation column: {col}")


expr_map = {}
for c in expression_df.columns:
    k = header_key(c)
    if k == "hugo_symbol":
        expr_map[c] = "hugo_symbol"
    elif k == "entrez_gene_id":
        expr_map[c] = "entrez_id"

expression_df = expression_df.rename(columns=expr_map)
#expression_df = coalesce_duplicate_columns(expression_df, "expression")

for col in ["hugo_symbol", "entrez_id"]:
    if col not in expression_df.columns:
        raise ValueError(f"Missing expected expression column: {col}")


# =========================
# CLEAN VALUES
# =========================

for df in [patient_df, sample_df, mutation_df, expression_df]:
    for c in df.columns:
        df[c] = df[c].map(norm_text)

patient_df["patient_id"] = patient_df["patient_id"].map(norm_key)
sample_df["patient_id"] = sample_df["patient_id"].map(norm_key)
sample_df["sample_id"] = sample_df["sample_id"].map(norm_key)
mutation_df["sample_id"] = mutation_df["sample_id"].map(norm_key)


# =========================
# EXPRESSION WIDE TO LONG
# =========================

# Only melt columns that look like TCGA sample IDs.
# Extra expression annotation columns such as gene_name, description, cytoband,
# gene_type, transcript_id, etc. will be ignored.
expr_sample_cols = [
    c for c in expression_df.columns
    if c not in ("hugo_symbol", "entrez_id")
    and is_tcga_sample_column(c)
]

ignored_expr_cols = [
    c for c in expression_df.columns
    if c not in ("hugo_symbol", "entrez_id") and c not in expr_sample_cols
]

if ignored_expr_cols:
    print("Ignoring non-sample expression columns:", ignored_expr_cols)

if not expr_sample_cols:
    raise ValueError(
        "No expression sample columns found. Expected columns that look like "
        "TCGA sample IDs, for example TCGA-XX-XXXX or TCGA_XX_XXXX."
    )

expression_long = expression_df.melt(
    id_vars=["hugo_symbol", "entrez_id"],
    value_vars=expr_sample_cols,
    var_name="sample_id",
    value_name="expression_value"
)

# Convert cleaned header sample IDs back to TCGA-style hyphen IDs so they can
# match sample_id values from the sample and mutation files.
expression_long["sample_id"] = expression_long["sample_id"].map(
    lambda x: norm_key(str(x).replace("_", "-")) if norm_key(x) else None
)
expression_long["hugo_symbol"] = expression_long["hugo_symbol"].map(norm_key)
expression_long["entrez_id"] = expression_long["entrez_id"].map(norm_key)
expression_long["expression_value"] = expression_long["expression_value"].map(norm_key)

expression_long = expression_long[
    expression_long["sample_id"].notna()
    & expression_long["hugo_symbol"].notna()
    & expression_long["entrez_id"].notna()
    & expression_long["expression_value"].notna()
].drop_duplicates(subset=["sample_id", "hugo_symbol", "entrez_id"])


# =========================
# CORE DATAFRAMES
# =========================

patient_core = patient_df[
    ["patient_id", "gender", "race", "ethnicity", "stage", "age", "cancer_type"]
].drop_duplicates()

sample_core = sample_df[
    ["sample_id", "patient_id", "tumor_tissue"]
].drop_duplicates()

mutation_core = mutation_df[
    [
        "sample_id", "hugo_symbol", "entrez_id", "chromosome",
        "hgvsp", "variant_classification", "codon_change",
        "start_pos", "end_pos", "dbsnp", "reference_allele",
        "tumor_allele1", "tumor_allele2", "polyphen_score", "variant_type"
    ]
].drop_duplicates()


# =========================
# CREATE ID MAPS
# =========================

race_vals = sorted({norm_key(v) for v in patient_core["race"] if norm_key(v)})
ethnicity_vals = sorted({norm_key(v) for v in patient_core["ethnicity"] if norm_key(v)})
disease_vals = sorted({norm_key(v) for v in patient_core["cancer_type"] if norm_key(v)})

stage_vals = sorted(
    {
        norm_key(r["stage"])
        for _, r in patient_core.iterrows()
        if norm_key(r["stage"])
    }
)

# Genes are unique by hugo_symbol + entrez_id.
# Chromosome is stored as an attribute, preferably from the mutation file.
gene_by_pair = {}

# 1. Add mutation genes first because mutation file usually has chromosome.
for _, r in mutation_core.iterrows():
    hugo = norm_gene_symbol(r["hugo_symbol"])
    entrez = norm_entrez_id(r["entrez_id"])
    chrom = norm_key(r["chromosome"])

    if hugo and entrez:
        pair = (hugo, entrez)

        # Prefer the first non-empty chromosome found from mutation data.
        if pair not in gene_by_pair:
            gene_by_pair[pair] = chrom
        elif gene_by_pair[pair] is None and chrom is not None:
            gene_by_pair[pair] = chrom

# 2. Add expression genes only if the hugo_symbol + entrez_id pair is new.
for _, r in expression_long.iterrows():
    hugo = norm_gene_symbol(r.get("hugo_symbol"))
    entrez = norm_entrez_id(r["entrez_id"])

    if hugo and entrez:
        pair = (hugo, entrez)

        if pair not in gene_by_pair:
            gene_by_pair[pair] = None

gene_vals = sorted(
    [
        (hugo, entrez, chrom)
        for (hugo, entrez), chrom in gene_by_pair.items()
    ],
    key=sort_key_tuple
)

mutation_vals = []

for _, r in mutation_core.iterrows():
    key = (
        norm_gene_symbol(r["hugo_symbol"]),
        norm_entrez_id(r["entrez_id"]),
        norm_key(r["chromosome"]),
        norm_key(r["hgvsp"]),
        norm_key(r["variant_classification"]),
        norm_key(r["codon_change"]),
        norm_key(r["start_pos"]),
        norm_key(r["end_pos"]),
        norm_key(r["dbsnp"]),
        norm_key(r["reference_allele"]),
        norm_key(r["tumor_allele1"]),
        norm_key(r["tumor_allele2"]),
        norm_key(r["polyphen_score"]),
        norm_key(r["variant_type"]),
    )

    mutation_vals.append(key)

mutation_vals = sorted(set(mutation_vals), key=sort_key_tuple)

race_id_map = make_id_map(race_vals)
ethnicity_id_map = make_id_map(ethnicity_vals)
disease_id_map = make_id_map(disease_vals)
stage_id_map = make_id_map(stage_vals)
gene_id_map = make_id_map(gene_vals)
mutation_id_map = make_id_map(mutation_vals)

# Expression and mutation rows should match genes by hugo_symbol + entrez_id.
gene_pair_id_map = {
    (hugo, entrez): gene_id
    for (hugo, entrez, chrom), gene_id in gene_id_map.items()
}


# =========================
# GENERATE SQL
# =========================

sql_parts = []

sql_parts.append(f"USE {DB_NAME};")
sql_parts.append("SET FOREIGN_KEY_CHECKS = 0;")
sql_parts.append("START TRANSACTION;")

sql_parts.append("""
-- Optional fresh reload
-- DELETE FROM expression;
-- DELETE FROM sample_mutation;
-- DELETE FROM mutation;
-- DELETE FROM sample;
-- DELETE FROM patient;
-- DELETE FROM gene;
-- DELETE FROM stage;
-- DELETE FROM disease;
-- DELETE FROM ethnicity;
-- DELETE FROM race;
""".strip())


# =========================
# LOOKUP TABLES WITH IDS
# =========================

race_rows = [
    [str(race_id_map[v]), sql_quote(v)]
    for v in race_vals
]

ethnicity_rows = [
    [str(ethnicity_id_map[v]), sql_quote(v)]
    for v in ethnicity_vals
]

disease_rows = [
    [str(disease_id_map[v]), sql_quote(v)]
    for v in disease_vals
]

stage_rows = [
    [str(stage_id_map[v]), sql_quote(v)]
    for v in stage_vals
]

gene_rows = [
    [str(gene_id_map[v]), sql_quote(v[0]), sql_quote(v[1]), sql_quote(v[2])]
    for v in gene_vals
]

sql_parts.extend(multirow_insert("race", ["race_id", "race"], race_rows))
sql_parts.extend(multirow_insert("ethnicity", ["ethnicity_id", "ethnicity"], ethnicity_rows))
sql_parts.extend(multirow_insert("disease", ["disease_id", "cancer_type"], disease_rows))
sql_parts.extend(multirow_insert("stage", ["stage_id", "stage"], stage_rows))
sql_parts.extend(multirow_insert("gene", ["gene_id", "hugo_symbol", "entrez_id", "chromosome"], gene_rows))


# =========================
# PATIENT
# =========================

patient_rows = []

for _, r in patient_core.iterrows():
    race_id = race_id_map.get(norm_key(r["race"]))
    ethnicity_id = ethnicity_id_map.get(norm_key(r["ethnicity"]))

    patient_rows.append([
        sql_quote(r["patient_id"]),
        sql_quote(r["gender"]),
        sql_number(race_id),
        sql_number(ethnicity_id)
    ])

sql_parts.extend(
    multirow_insert(
        "patient",
        ["patient_id", "gender", "race_id", "ethnicity_id"],
        patient_rows
    )
)


# =========================
# SAMPLE
# =========================

sample_merged = sample_core.merge(
    patient_core[["patient_id", "cancer_type", "age", "stage"]],
    on="patient_id",
    how="left"
)

sample_rows = []

for _, r in sample_merged.iterrows():
    disease_id = disease_id_map.get(norm_key(r["cancer_type"]))
    stage_key = (norm_key(r["stage"]))
    stage_id = stage_id_map.get(stage_key)

    sample_rows.append([
        sql_quote(r["sample_id"]),
        sql_quote(r["patient_id"]),
        sql_number(disease_id),
        sql_number(r["age"]),
        sql_number(stage_id),
        sql_quote(r["tumor_tissue"])
    ])

sql_parts.extend(
    multirow_insert(
        "sample",
        ["sample_id", "patient_id", "disease_id", "age", "stage_id", "tumor_tissue"],
        sample_rows
    )
)


# =========================
# MUTATION
# =========================

mutation_rows = []

for key in mutation_vals:
    (
        hugo_symbol,
        entrez_id,
        chromosome,
        hgvsp,
        variant_classification,
        codon_change,
        start_pos,
        end_pos,
        dbsnp,
        reference_allele,
        tumor_allele1,
        tumor_allele2,
        polyphen_score,
        variant_type,
    ) = key

    gene_id = gene_pair_id_map.get((hugo_symbol, entrez_id))
    mutation_id = mutation_id_map[key]

    mutation_rows.append([
        str(mutation_id),
        sql_number(gene_id),
        sql_quote(hgvsp),
        sql_quote(variant_classification),
        sql_quote(codon_change),
        sql_number(start_pos),
        sql_number(end_pos),
        sql_quote(dbsnp),
        sql_quote(reference_allele),
        sql_quote(tumor_allele1),
        sql_quote(tumor_allele2),
        sql_quote(polyphen_score),
        sql_quote(variant_type)
    ])

sql_parts.extend(
    multirow_insert(
        "mutation",
        [
            "mutation_id",
            "gene_id",
            "HGVSp",
            "variantClassification",
            "codon_change",
            "start_pos",
            "end_pos",
            "dbSNP",
            "reference_allele",
            "tumor_allele1",
            "tumor_allele2",
            "polyphen_score",
            "varian_type"
        ],
        mutation_rows
    )
)


# =========================
# SAMPLE_MUTATION
# =========================

sample_mutation_pairs = set()

for _, r in mutation_core.iterrows():
    mutation_key = (
        norm_gene_symbol(r["hugo_symbol"]),
        norm_entrez_id(r["entrez_id"]),
        norm_key(r["chromosome"]),
        norm_key(r["hgvsp"]),
        norm_key(r["variant_classification"]),
        norm_key(r["codon_change"]),
        norm_key(r["start_pos"]),
        norm_key(r["end_pos"]),
        norm_key(r["dbsnp"]),
        norm_key(r["reference_allele"]),
        norm_key(r["tumor_allele1"]),
        norm_key(r["tumor_allele2"]),
        norm_key(r["polyphen_score"]),
        norm_key(r["variant_type"]),
    )

    mutation_id = mutation_id_map.get(mutation_key)
    sample_id = norm_key(r["sample_id"])

    if sample_id and mutation_id:
        sample_mutation_pairs.add((sample_id, mutation_id))

sample_mutation_list = sorted(sample_mutation_pairs)

sample_mutation_id_map = {
    pair: i for i, pair in enumerate(sample_mutation_list, start=1)
}

sample_mutation_rows = [
    [
        str(sample_mutation_id_map[(sample_id, mutation_id)]),
        sql_quote(sample_id),
        str(mutation_id)
    ]
    for sample_id, mutation_id in sample_mutation_list
]

sql_parts.extend(
    multirow_insert(
        "sample_mutation",
        ["sample_mutation_id", "sample_id", "mutation_id"],
        sample_mutation_rows
    )
)


# =========================
# EXPRESSION
# =========================

expression_rows = []

for _, r in expression_long.iterrows():
    sample_id = norm_key(r["sample_id"])
    hugo_symbol = norm_gene_symbol(r["hugo_symbol"])
    entrez_id = norm_entrez_id(r["entrez_id"])

    gene_id = gene_pair_id_map.get((hugo_symbol, entrez_id))

    if sample_id and gene_id:
        expression_rows.append([
            sql_quote(sample_id),
            str(gene_id),
            sql_number(r["expression_value"])
        ])

sql_parts.extend(
    multirow_insert(
        "expression",
        ["sample_id", "gene_id", "expressionValue"],
        expression_rows
    )
)


# =========================
# RESET AUTO_INCREMENT
# =========================

sql_parts.append(f"ALTER TABLE race AUTO_INCREMENT = {len(race_vals) + 1};")
sql_parts.append(f"ALTER TABLE ethnicity AUTO_INCREMENT = {len(ethnicity_vals) + 1};")
sql_parts.append(f"ALTER TABLE disease AUTO_INCREMENT = {len(disease_vals) + 1};")
sql_parts.append(f"ALTER TABLE stage AUTO_INCREMENT = {len(stage_vals) + 1};")
sql_parts.append(f"ALTER TABLE gene AUTO_INCREMENT = {len(gene_vals) + 1};")
sql_parts.append(f"ALTER TABLE mutation AUTO_INCREMENT = {len(mutation_vals) + 1};")

sql_parts.append("COMMIT;")
sql_parts.append("SET FOREIGN_KEY_CHECKS = 1;")


# =========================
# WRITE SQL FILE
# =========================

Path(OUTPUT_SQL).write_text("\n\n".join(sql_parts), encoding="utf-8")

print(f"SQL written to {OUTPUT_SQL}")
print(f"Race rows: {len(race_vals)}")
print(f"Ethnicity rows: {len(ethnicity_vals)}")
print(f"Disease rows: {len(disease_vals)}")
print(f"Stage rows: {len(stage_vals)}")
print(f"Gene rows: {len(gene_vals)}")
print(f"Unique gene pairs: {len(gene_pair_id_map)}")
print(f"Mutation rows: {len(mutation_vals)}")
print(f"Sample-mutation rows: {len(sample_mutation_rows)}")
print(f"Expression rows: {len(expression_rows)}")
