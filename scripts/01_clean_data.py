import pandas as pd

# -----------------------------
# File paths
# -----------------------------
input_file = "data_mrna_seq_v2_rsem.txt"
output_file = "expression_data_cleaned1.tsv"

# -----------------------------
# 1. Read the tab-delimited file
# -----------------------------
# Skip metadata/comment rows starting with '#'
# Treat "." as a missing value
df = pd.read_csv(
    input_file,
    sep="\t",
    comment="#",
    dtype=str,
    keep_default_na=True,
    na_values=["", " ", ".", "NA", "N/A", "na", "null", "NULL", "None"]
)

# -----------------------------
# 2. Preserve column-name capitalization
# -----------------------------
# This only removes leading/trailing spaces from column names.
# It does NOT convert names to lowercase.
# It does NOT replace special characters.
df.columns = df.columns.str.strip()

# -----------------------------
# 3. Trim whitespace in all string cells
# -----------------------------
df = df.apply(lambda col: col.str.strip() if col.dtype == "object" else col)

# -----------------------------
# 4. Standardize null-like values again after trimming
# -----------------------------
# "." is included here too, in case it had surrounding spaces before trimming
null_like_values = {
    "": pd.NA,
    " ": pd.NA,
    ".": pd.NA,
    "na": pd.NA,
    "n/a": pd.NA,
    "null": pd.NA,
    "none": pd.NA,
    "unknown": pd.NA
}

df = df.replace(null_like_values)

# -----------------------------
# 5. Drop columns that have a header but no values
# -----------------------------
# This removes columns where every cell is missing/blank after cleanup.
columns_before = len(df.columns)
empty_columns = df.columns[df.isna().all()].tolist()
df = df.dropna(axis=1, how="all")
columns_after = len(df.columns)
empty_columns_removed = columns_before - columns_after

# -----------------------------
# 6. Remove fully duplicate rows
# -----------------------------
rows_before = len(df)
df = df.drop_duplicates()
rows_after = len(df)
duplicates_removed = rows_before - rows_after

# -----------------------------
# 7. Optional: remove duplicate patient_id rows
#    Keep the first occurrence only
# -----------------------------
patient_id_duplicates_removed = 0

if "patient_id" in df.columns:
    before_pid = len(df)
    df = df.drop_duplicates(subset=["patient_id"], keep="first")
    after_pid = len(df)
    patient_id_duplicates_removed = before_pid - after_pid

# -----------------------------
# 8. Save cleaned file as tab-delimited
# -----------------------------
# na_rep="" writes nulls as blank cells
df.to_csv(output_file, sep="\t", index=False, na_rep="")

# -----------------------------
# 9. Print summary
# -----------------------------
print("Cleaning complete.")
print(f"Input file:  {input_file}")
print(f"Output file: {output_file}")
print(f"Rows after cleaning: {len(df)}")
print(f"Empty columns removed: {empty_columns_removed}")
print(f"Full duplicate rows removed: {duplicates_removed}")
print(f"Duplicate patient_id rows removed: {patient_id_duplicates_removed}")

if empty_columns:
    print("\nColumns removed because they were completely empty:")
    for col in empty_columns:
        print(f"- {col}")

print("\nNull values by column:")
print(df.isna().sum())