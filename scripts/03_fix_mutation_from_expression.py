#!/usr/bin/env python3
"""
Fix mutation-file Entrez IDs using HUGO symbols from an expression file.

This version uses hardcoded file paths instead of command-line arguments.
Edit the paths in the CONFIG section before running.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import pandas as pd


# =========================
# CONFIG: EDIT THESE PATHS
# =========================
MUTATION_FILE = r"C:\Users\Hangt\Desktop\Bioinformatics\database\cancer database project\cleaned_data\mutation_mapped1.tsv"
EXPRESSION_FILE = r"C:\Users\Hangt\Desktop\Bioinformatics\database\cancer database project\cleaned_data\expression_mapped1.tsv"
OUTPUT_FILE = r"mutation_mapped_fixed2.tsv"

# Optional audit files. Set to None if you do not want them.
AUDIT_OUTPUT_FILE = r"C:\Users\Hangt\Desktop\Bioinformatics\database\cancer database project\changed_rows.csv"
MAPPING_AUDIT_OUTPUT_FILE = r"C:\Users\Hangt\Desktop\Bioinformatics\database\cancer database project\mapping_audit.csv"

# Optional Excel sheet names. Leave as None for CSV/TSV files.
MUTATION_SHEET = None
EXPRESSION_SHEET = None

# False = only replace missing/0 Entrez IDs in mutation file
# True = also replace existing non-zero Entrez IDs when a mapping exists
REPLACE_NONZERO = False


def normalize_symbol(value) -> Optional[str]:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "null"}:
        return None
    return text


def normalize_entrez(value) -> Optional[int]:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "null"}:
        return None
    try:
        numeric = int(float(text))
    except ValueError:
        return None
    if numeric == 0:
        return None
    return numeric


def infer_separator(path: Path) -> Optional[str]:
    suffix = path.suffix.lower()
    if suffix == ".tsv":
        return "\t"
    if suffix == ".csv":
        return ","
    return None


def read_table(path: Path, sheet_name: Optional[str] = None) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name=sheet_name)
    sep = infer_separator(path)
    if sep is None:
        return pd.read_csv(path, sep=None, engine="python")
    return pd.read_csv(path, sep=sep)


def write_table(df: pd.DataFrame, path: Path) -> None:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        df.to_csv(path, index=False)
    elif suffix == ".tsv":
        df.to_csv(path, sep="\t", index=False)
    elif suffix in {".xlsx", ".xls"}:
        df.to_excel(path, index=False)
    else:
        raise ValueError(f"Unsupported output format: {path.suffix}")


def find_column(columns: Iterable[str], candidates: Iterable[str]) -> str:
    lowered = {str(col).strip().lower(): col for col in columns}
    for candidate in candidates:
        key = candidate.strip().lower()
        if key in lowered:
            return lowered[key]
    raise KeyError(
        f"Could not find any of these columns: {list(candidates)}. "
        f"Available columns: {list(columns)}"
    )


def build_expression_mapping(
    expr_df: pd.DataFrame,
    expr_symbol_col: str,
    expr_entrez_col: str,
) -> Tuple[Dict[str, int], pd.DataFrame]:
    tmp = expr_df[[expr_symbol_col, expr_entrez_col]].copy()
    tmp["symbol_norm"] = tmp[expr_symbol_col].map(normalize_symbol)
    tmp["entrez_norm"] = tmp[expr_entrez_col].map(normalize_entrez)
    tmp = tmp.dropna(subset=["symbol_norm", "entrez_norm"])

    grouped = (
        tmp.groupby("symbol_norm")["entrez_norm"]
        .agg(lambda s: sorted(set(int(x) for x in s if pd.notna(x))))
        .reset_index()
        .rename(columns={"entrez_norm": "unique_entrez_ids"})
    )
    grouped["n_unique_entrez_ids"] = grouped["unique_entrez_ids"].map(len)
    grouped["status"] = grouped["n_unique_entrez_ids"].map(
        lambda n: "usable" if n == 1 else "ambiguous"
    )

    usable = grouped[grouped["status"] == "usable"].copy()
    mapping = {
        row["symbol_norm"]: int(row["unique_entrez_ids"][0])
        for _, row in usable.iterrows()
    }
    return mapping, grouped


def update_mutation_entrez_ids(
    mut_df: pd.DataFrame,
    mut_symbol_col: str,
    mut_entrez_col: str,
    mapping: Dict[str, int],
    replace_nonzero: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    out = mut_df.copy()
    out["_symbol_norm"] = out[mut_symbol_col].map(normalize_symbol)
    out["_entrez_norm"] = out[mut_entrez_col].map(normalize_entrez)
    out["_mapped_entrez"] = out["_symbol_norm"].map(mapping)

    if replace_nonzero:
        to_update = out["_mapped_entrez"].notna()
    else:
        to_update = out["_entrez_norm"].isna() & out["_mapped_entrez"].notna()

    audit = out.loc[
        to_update,
        [mut_symbol_col, mut_entrez_col, "_mapped_entrez"]
    ].copy()
    audit = audit.rename(
        columns={
            mut_symbol_col: "hugo_symbol",
            mut_entrez_col: "old_entrez_id",
            "_mapped_entrez": "new_entrez_id",
        }
    )

    out.loc[to_update, mut_entrez_col] = out.loc[to_update, "_mapped_entrez"].astype("Int64")
    out = out.drop(columns=["_symbol_norm", "_entrez_norm", "_mapped_entrez"])
    return out, audit


def main() -> None:
    mutation_path = Path(MUTATION_FILE)
    expression_path = Path(EXPRESSION_FILE)
    output_path = Path(OUTPUT_FILE)
    audit_output_path = Path(AUDIT_OUTPUT_FILE) if AUDIT_OUTPUT_FILE else None
    mapping_audit_output_path = Path(MAPPING_AUDIT_OUTPUT_FILE) if MAPPING_AUDIT_OUTPUT_FILE else None

    mut_df = read_table(mutation_path, sheet_name=MUTATION_SHEET)
    expr_df = read_table(expression_path, sheet_name=EXPRESSION_SHEET)

    mut_symbol_col = find_column(mut_df.columns, ["hugo_symbol", "Hugo_Symbol"])
    mut_entrez_col = find_column(mut_df.columns, ["entrez_gene_id", "Entrez_Gene_Id", "entrez_id"])
    expr_symbol_col = find_column(expr_df.columns, ["Hugo_Symbol", "hugo_symbol"])
    expr_entrez_col = find_column(expr_df.columns, ["Entrez_Gene_Id", "entrez_gene_id", "entrez_id"])

    mapping, mapping_audit = build_expression_mapping(
        expr_df=expr_df,
        expr_symbol_col=expr_symbol_col,
        expr_entrez_col=expr_entrez_col,
    )

    fixed_mut_df, row_audit = update_mutation_entrez_ids(
        mut_df=mut_df,
        mut_symbol_col=mut_symbol_col,
        mut_entrez_col=mut_entrez_col,
        mapping=mapping,
        replace_nonzero=REPLACE_NONZERO,
    )

    write_table(fixed_mut_df, output_path)

    if audit_output_path:
        write_table(row_audit, audit_output_path)

    if mapping_audit_output_path:
        write_table(mapping_audit, mapping_audit_output_path)

    summary = {
        "mutation_rows": int(len(mut_df)),
        "expression_rows": int(len(expr_df)),
        "usable_symbol_to_entrez_mappings": int((mapping_audit["status"] == "usable").sum()),
        "ambiguous_symbols_in_expression": int((mapping_audit["status"] == "ambiguous").sum()),
        "mutation_rows_updated": int(len(row_audit)),
        "replace_nonzero_enabled": bool(REPLACE_NONZERO),
        "mutation_symbol_column": mut_symbol_col,
        "mutation_entrez_column": mut_entrez_col,
        "expression_symbol_column": expr_symbol_col,
        "expression_entrez_column": expr_entrez_col,
        "output_file": str(output_path),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
