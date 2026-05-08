#!/usr/bin/env python3

import csv
from collections import defaultdict
from copy import deepcopy


# =========================================================
# SETTINGS: EDIT THESE VALUES IN PYCHARM, THEN CLICK RUN
# =========================================================

HGNC_FILE = r"hgnc_complete_set.txt"
MUTATION_FILE = r"mutations_data_cleaned.tsv"
EXPRESSION_FILE = r"expression_data_cleaned1.tsv"

MUTATION_OUTPUT = r"mutation_mapped1.tsv"
EXPRESSION_OUTPUT = r"expression_mapped1.tsv"

# Column names in your mutation file
MUTATION_SYMBOL_COL = "Hugo_Symbol"
MUTATION_ENTREZ_COL = "Entrez_Gene_Id"

# Column names in your expression file
EXPRESSION_SYMBOL_COL = "Hugo_Symbol"
EXPRESSION_ENTREZ_COL = "Entrez_Gene_Id"

# If True, replace the original symbol / entrez columns
# If False, keep originals and add mapped_* audit columns
OVERWRITE_ORIGINAL = True


# =========================================================
# HELPER FUNCTIONS
# =========================================================

def normalize_value(x):
    if x is None:
        return ""
    x = str(x).strip()
    if x in {"", "nan", "None", "NULL", "NA", "N/A"}:
        return ""
    return x


def normalize_symbol(x):
    return normalize_value(x).upper()


def normalize_entrez(x):
    x = normalize_value(x)
    if not x:
        return ""
    return x


def split_multi_value_field(value):
    value = normalize_value(value)
    if not value:
        return []

    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        value = value[1:-1]

    parts = [p.strip() for p in value.split("|")]
    return [p for p in parts if p]


def read_tsv(path):
    with open(path, "r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    return fieldnames, rows


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_hgnc_maps(hgnc_rows):
    approved_rows = [r for r in hgnc_rows if normalize_value(r.get("status")) == "Approved"]

    entrez_map = {}
    symbol_map = {}
    prev_symbol_map = defaultdict(list)
    alias_symbol_map = defaultdict(list)

    for row in approved_rows:
        approved_symbol = normalize_symbol(row.get("symbol"))
        entrez_id = normalize_entrez(row.get("entrez_id"))

        if approved_symbol:
            symbol_map[approved_symbol] = row

        if entrez_id and entrez_id not in entrez_map:
            entrez_map[entrez_id] = row

        for prev in split_multi_value_field(row.get("prev_symbol")):
            prev_norm = normalize_symbol(prev)
            if prev_norm:
                prev_symbol_map[prev_norm].append(row)

        for alias in split_multi_value_field(row.get("alias_symbol")):
            alias_norm = normalize_symbol(alias)
            if alias_norm:
                alias_symbol_map[alias_norm].append(row)

    return approved_rows, entrez_map, symbol_map, prev_symbol_map, alias_symbol_map


def choose_unique_match(matches):
    if not matches:
        return None

    unique = {}
    for row in matches:
        sym = normalize_symbol(row.get("symbol"))
        if sym:
            unique[sym] = row

    if len(unique) == 1:
        return next(iter(unique.values()))
    return None


def map_row(
    row,
    symbol_col,
    entrez_col,
    entrez_map,
    symbol_map,
    prev_symbol_map,
    alias_symbol_map
):
    original_symbol = normalize_value(row.get(symbol_col))
    original_entrez = normalize_value(row.get(entrez_col))

    symbol_norm = normalize_symbol(original_symbol)
    entrez_norm = normalize_entrez(original_entrez)

    result = {
        "original_symbol": original_symbol,
        "original_entrez_id": original_entrez,
        "mapped_symbol": original_symbol,
        "mapped_entrez_id": original_entrez,
        "mapped_hgnc_id": "",
        "mapping_source": "",
        "mapping_status": "",
    }

    entrez_hit = entrez_map.get(entrez_norm) if entrez_norm else None
    symbol_hit = symbol_map.get(symbol_norm) if symbol_norm else None
    prev_candidates = prev_symbol_map.get(symbol_norm, []) if symbol_norm else []
    alias_candidates = alias_symbol_map.get(symbol_norm, []) if symbol_norm else []
    prev_hit = choose_unique_match(prev_candidates) if symbol_norm else None
    alias_hit = choose_unique_match(alias_candidates) if symbol_norm else None

    if entrez_hit is not None:
        symbol_based_hit = symbol_hit or prev_hit or alias_hit
        if symbol_based_hit is not None:
            entrez_sym = normalize_symbol(entrez_hit.get("symbol"))
            symbol_sym = normalize_symbol(symbol_based_hit.get("symbol"))
            if entrez_sym != symbol_sym:
                result["mapping_status"] = "conflict"
                result["mapping_source"] = "entrez_vs_symbol"
                return result

    chosen = None
    source = ""

    if entrez_hit is not None:
        chosen = entrez_hit
        source = "entrez"
    elif symbol_hit is not None:
        chosen = symbol_hit
        source = "symbol"
    elif prev_hit is not None:
        chosen = prev_hit
        source = "prev_symbol"
    elif prev_candidates:
        result["mapping_status"] = "ambiguous"
        result["mapping_source"] = "prev_symbol"
        return result
    elif alias_hit is not None:
        chosen = alias_hit
        source = "alias_symbol"
    elif alias_candidates:
        result["mapping_status"] = "ambiguous"
        result["mapping_source"] = "alias_symbol"
        return result
    else:
        result["mapping_status"] = "missing"
        result["mapping_source"] = "none"
        return result

    approved_symbol = normalize_value(chosen.get("symbol"))
    approved_entrez = normalize_value(chosen.get("entrez_id"))
    hgnc_id = normalize_value(chosen.get("hgnc_id"))

    result["mapped_symbol"] = approved_symbol or original_symbol
    result["mapped_entrez_id"] = approved_entrez or original_entrez
    result["mapped_hgnc_id"] = hgnc_id
    result["mapping_source"] = source

    if (
        normalize_symbol(original_symbol) == normalize_symbol(approved_symbol)
        and normalize_entrez(original_entrez) == normalize_entrez(approved_entrez)
    ):
        result["mapping_status"] = "exact"
    else:
        result["mapping_status"] = "updated"

    return result


def process_gene_file(
    input_path,
    output_path,
    hgnc_rows,
    symbol_col,
    entrez_col,
    overwrite_original=False
):
    fieldnames, rows = read_tsv(input_path)

    if symbol_col not in fieldnames:
        raise ValueError(f"Symbol column '{symbol_col}' not found in file: {input_path}")
    if entrez_col not in fieldnames:
        raise ValueError(f"Entrez column '{entrez_col}' not found in file: {input_path}")

    _, entrez_map, symbol_map, prev_symbol_map, alias_symbol_map = build_hgnc_maps(hgnc_rows)

    audit_cols = [
        "original_symbol",
        "original_entrez_id",
        "mapped_symbol",
        "mapped_entrez_id",
        "mapped_hgnc_id",
        "mapping_source",
        "mapping_status",
    ]

    output_rows = []
    for row in rows:
        out = deepcopy(row)

        mapped = map_row(
            row=out,
            symbol_col=symbol_col,
            entrez_col=entrez_col,
            entrez_map=entrez_map,
            symbol_map=symbol_map,
            prev_symbol_map=prev_symbol_map,
            alias_symbol_map=alias_symbol_map,
        )

        # If a gene symbol is present/mapped but Entrez could not be mapped,
        # some input files use 0 as a placeholder. Write that as blank instead.
        # The original_entrez_id audit column still preserves the original value.
        if normalize_value(mapped.get("mapped_symbol")) and normalize_entrez(mapped.get("mapped_entrez_id")) == "0":
            mapped["mapped_entrez_id"] = ""

        for col in audit_cols:
            out[col] = mapped[col]

        if overwrite_original:
            out[symbol_col] = mapped["mapped_symbol"]
            out[entrez_col] = mapped["mapped_entrez_id"]

        output_rows.append(out)

    final_fieldnames = fieldnames[:]
    for col in audit_cols:
        if col not in final_fieldnames:
            final_fieldnames.append(col)

    write_tsv(output_path, final_fieldnames, output_rows)


def main():
    print("Reading HGNC reference file...")
    _, hgnc_rows = read_tsv(HGNC_FILE)

    print("Processing mutation file...")
    process_gene_file(
        input_path=MUTATION_FILE,
        output_path=MUTATION_OUTPUT,
        hgnc_rows=hgnc_rows,
        symbol_col=MUTATION_SYMBOL_COL,
        entrez_col=MUTATION_ENTREZ_COL,
        overwrite_original=OVERWRITE_ORIGINAL,
    )

    print("Processing expression file...")
    process_gene_file(
        input_path=EXPRESSION_FILE,
        output_path=EXPRESSION_OUTPUT,
        hgnc_rows=hgnc_rows,
        symbol_col=EXPRESSION_SYMBOL_COL,
        entrez_col=EXPRESSION_ENTREZ_COL,
        overwrite_original=OVERWRITE_ORIGINAL,
    )

    print("Done.")
    print(f"Mutation output:   {MUTATION_OUTPUT}")
    print(f"Expression output: {EXPRESSION_OUTPUT}")


if __name__ == "__main__":
    main()