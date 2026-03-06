#!/usr/bin/env python3
"""
Clean and organize clinical data for BRCA piRNA pipeline.

Key discovery: The clinical files are CROSS-MATCHED with expression files:
  - yyfbatch2_clinical.csv RNA_ID2 -> yyfbatch1_processed.csv row names (exact match)
  - yyfbatch1_clinical.csv RNA_ID  -> yyfbatch2_processed.csv row names (prefix match)

This script:
  1. Builds ID mapping between clinical and expression data
  2. Standardizes staging columns (AJCC Overall, T, N, M) for both clinical and pathologic
  3. Derives molecular subtype from ER/PR/HER2
  4. Outputs clean CSVs where SampleID = expression row name
"""

import csv
import re
import os
import sys

# ---------------------------------------------------------------------------
# 0. Helpers
# ---------------------------------------------------------------------------

def read_csv(path):
    """Read CSV handling encoding issues."""
    for enc in ("utf-8", "latin-1", "cp1252"):
        try:
            with open(path, "r", encoding=enc) as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            return rows
        except UnicodeDecodeError:
            continue
    raise RuntimeError(f"Cannot decode {path}")


def get_expression_ids(path):
    """Get row names (first column) from expression CSV."""
    ids = []
    for enc in ("utf-8", "latin-1", "cp1252"):
        try:
            with open(path, "r", encoding=enc) as f:
                reader = csv.reader(f)
                header = next(reader)
                for row in reader:
                    ids.append(row[0])
            return ids
        except UnicodeDecodeError:
            continue
    raise RuntimeError(f"Cannot decode {path}")


def clean_str(s):
    """Strip whitespace and normalize NA."""
    if s is None:
        return ""
    s = s.strip()
    if s.upper() in ("NA", "N/A", "NONE", "", "NOT APPLICABLE", "PM: NOT APPLICABLE"):
        return ""
    return s


# ---------------------------------------------------------------------------
# 1. Standardize staging
# ---------------------------------------------------------------------------

def parse_overall_stage(raw):
    """Parse Best CS/AJCC Stage into standard stage group.

    Returns (stage_standard, stage_group) e.g. ('Stage IIB', 'Stage II')
    """
    s = clean_str(raw).upper().replace(" ", "")

    if not s or s in ("NA",):
        return "", ""

    # Strip "STAGE" prefix if present
    s = s.replace("STAGE", "").strip()

    # Map numeric-letter combos to roman numeral stages
    stage_map = {
        "0": ("Stage 0", "Stage 0"),
        "1": ("Stage I", "Stage I"),
        "1A": ("Stage IA", "Stage I"),
        "1B": ("Stage IB", "Stage I"),
        "1C": ("Stage IC", "Stage I"),
        "2A": ("Stage IIA", "Stage II"),
        "2B": ("Stage IIB", "Stage II"),
        "3A": ("Stage IIIA", "Stage III"),
        "3B": ("Stage IIIB", "Stage III"),
        "3C": ("Stage IIIC", "Stage III"),
        "4": ("Stage IV", "Stage IV"),
        "4A": ("Stage IVA", "Stage IV"),
        "4B": ("Stage IVB", "Stage IV"),
    }

    if s in stage_map:
        return stage_map[s]

    # Already roman numeral format
    roman_map = {
        "0": "Stage 0",
        "I": "Stage I", "IA": "Stage IA", "IB": "Stage IB", "IC": "Stage IC",
        "II": "Stage II", "IIA": "Stage IIA", "IIB": "Stage IIB",
        "III": "Stage III", "IIIA": "Stage IIIA", "IIIB": "Stage IIIB", "IIIC": "Stage IIIC",
        "IV": "Stage IV", "IVA": "Stage IVA", "IVB": "Stage IVB",
    }
    if s in roman_map:
        full = roman_map[s]
        group = re.sub(r"[ABC]$", "", full)
        return full, group

    return raw.strip(), ""


def parse_tnm(raw):
    """Parse a T, N, or M staging value into a standardized category.

    Handles formats like: cT2, pT1c, c2, p0I-, pIS, pTis(DCIS), pT1b(m), etc.
    Returns (value, prefix) where prefix is 'c' (clinical) or 'p' (pathologic).
    """
    s = clean_str(raw).strip()
    if not s:
        return "", ""

    # Determine prefix (clinical vs pathologic)
    # Handle yp (post-neoadjuvant pathologic), yc, etc.
    prefix = ""
    val = s
    val_lower = val.lower()
    if val_lower.startswith("yp"):
        prefix = "p"
        val = val[2:]
    elif val_lower.startswith("yc"):
        prefix = "c"
        val = val[2:]
    elif val_lower.startswith("c"):
        prefix = "c"
        val = val[1:]
    elif val_lower.startswith("p"):
        prefix = "p"
        val = val[1:]

    # Remove leading T/N/M letter if present (e.g., cT2 -> 2, pN0 -> 0)
    val_upper = val.upper()
    if val_upper.startswith(("T", "N", "M")):
        val = val[1:]

    # Clean up: remove parenthetical notes, trailing whitespace
    val = re.sub(r"\s+", "", val)
    # Keep the base value but remove annotations like (i+), (sn), (m), (DCIS)
    val_clean = re.sub(r"\(.*?\)", "", val)

    return val_clean.strip(), prefix


def standardize_t_stage(raw_clinical, raw_pathologic):
    """Standardize T stage. Prefer pathologic if available, else clinical.

    Returns (T_value, T_category) where T_category is for grouping:
      Tis, T0, T1, T2, T3, T4
    """
    # Try pathologic first, then clinical
    for raw in [raw_pathologic, raw_clinical]:
        val, _ = parse_tnm(raw)
        if val:
            val_upper = val.upper()
            # Categorize
            if "IS" in val_upper or "DCIS" in val_upper:
                return f"T{val}", "Tis"
            elif val_upper.startswith("X") or val_upper == "X":
                return "TX", "TX"
            else:
                # Extract the leading digit
                m = re.match(r"(\d)", val)
                if m:
                    digit = m.group(1)
                    return f"T{val}", f"T{digit}"
            return f"T{val}", f"T{val}"
    return "", ""


def standardize_n_stage(raw_clinical, raw_pathologic):
    """Standardize N stage. Returns (N_value, N_category)."""
    for raw in [raw_pathologic, raw_clinical]:
        val, _ = parse_tnm(raw)
        if val:
            val_upper = val.upper()
            if val_upper.startswith("X") or val_upper == "X":
                return "NX", "NX"
            m = re.match(r"(\d)", val)
            if m:
                digit = m.group(1)
                return f"N{val}", f"N{digit}"
            return f"N{val}", f"N{val}"
    return "", ""


def standardize_m_stage(raw_clinical, raw_pathologic):
    """Standardize M stage. Returns (M_value, M_category)."""
    for raw in [raw_pathologic, raw_clinical]:
        val, _ = parse_tnm(raw)
        if val:
            val_upper = val.upper()
            if val_upper.startswith("X") or val_upper == "X":
                return "MX", "MX"
            m = re.match(r"(\d)", val)
            if m:
                digit = m.group(1)
                return f"M{val}", f"M{digit}"
            return f"M{val}", f"M{val}"
    return "", ""


# ---------------------------------------------------------------------------
# 2. Derive molecular subtype from ER/PR/HER2
# ---------------------------------------------------------------------------

def derive_subtype(er_raw, pr_raw, her2_raw):
    """Derive molecular subtype from receptor status.

    Returns (Subtype, ER_status, PR_status, HER2_status)
    """
    def parse_receptor(s):
        s = clean_str(s).lower()
        if not s:
            return "unknown"
        # Handle detailed formats: "ER positive (20%)", "Positive (>95%)",
        # "Hercept = Equivocal (2+) ; FISH = Not amplified", "Negative (1+)"
        # Check for FISH amplification override
        if "fish" in s:
            if "not amplified" in s or "non-amplified" in s:
                return "negative"
            elif "amplified" in s:
                return "positive"
        if any(w in s for w in ("positive", "postive", "pos ")):
            return "positive"
        if "negative" in s or s == "neg" or s == "-":
            return "negative"
        if "equivocal" in s or "2+" in s:
            return "equivocal"
        if s in ("+",):
            return "positive"
        return "unknown"

    er = parse_receptor(er_raw)
    pr = parse_receptor(pr_raw)
    her2 = parse_receptor(her2_raw)

    # Subtype classification
    if er == "positive" or pr == "positive":
        if her2 == "positive":
            subtype = "Luminal B"
        else:
            # Luminal A vs B ideally needs Ki-67, default to Luminal A
            subtype = "Luminal A"
    elif her2 == "positive":
        subtype = "HER2+"
    elif er == "negative" and pr == "negative" and her2 == "negative":
        subtype = "Triple-negative"
    else:
        subtype = "Unclassified"

    return subtype, er, pr, her2


def parse_alive_status(raw):
    """Parse Alive Yes/No into standardized OS_status."""
    s = clean_str(raw).upper()
    if s in ("YES", "Y"):
        return 0  # alive = censored
    elif s in ("NO", "N"):
        return 1  # dead = event
    return ""  # unknown


def parse_os_months(raw):
    """Parse OS (Mo) into numeric."""
    s = clean_str(raw)
    try:
        return round(float(s), 2)
    except (ValueError, TypeError):
        return ""


# ---------------------------------------------------------------------------
# 3. Build ID mappings
# ---------------------------------------------------------------------------

def build_batch1expr_to_batch2clin(batch2_clin):
    """batch1 expression row names (vCa1_merged) <- batch2 clinical RNA_ID2."""
    mapping = {}
    for row in batch2_clin:
        rna_id2 = clean_str(row.get("RNA_ID2", ""))
        if rna_id2:
            mapping[rna_id2] = row
    return mapping


def build_batch2expr_to_batch1clin(batch1_clin, batch2_expr_ids):
    """batch2 expression row names (LA1_sample1_R1) <- batch1 clinical RNA_ID.

    Mapping logic:
      - expression prefix before _sample or _bksample -> RNA_ID for most
      - B{n}_sample* -> Be{n} in clinical
      - b{digits}A_bksample* -> source_ID M2200{digits}A -> BB sample
      - c{digits}A_bksample* -> source_ID M2200{digits}A -> BC sample
    """
    # Build lookup by RNA_ID
    by_rna_id = {}
    for row in batch1_clin:
        rid = clean_str(row.get("RNA_ID", ""))
        if rid:
            by_rna_id[rid] = row

    # Build lookup by source_ID (for BioBank b*/c* samples)
    by_source_id = {}
    for row in batch1_clin:
        sid = clean_str(row.get("source_ID", ""))
        if sid:
            by_source_id[sid] = row

    mapping = {}
    unmatched = []
    for expr_id in batch2_expr_ids:
        # Extract prefix before _sample or _bksample
        prefix = re.sub(r"_sample.*|_bksample.*", "", expr_id)

        # Direct RNA_ID match (LA1, HR1, TN1, LB1, He1, etc.)
        if prefix in by_rna_id:
            mapping[expr_id] = by_rna_id[prefix]
            continue

        # B{n} -> Be{n} mapping for benign
        m_benign = re.match(r"^B(\d+)$", prefix)
        if m_benign:
            be_id = f"Be{m_benign.group(1)}"
            if be_id in by_rna_id:
                mapping[expr_id] = by_rna_id[be_id]
                continue

        # b{digits}A -> M2200{digits}A (BioBank benign -> BB samples)
        m_bk_benign = re.match(r"^b(\d+)A$", prefix)
        if m_bk_benign:
            source_id = f"M2200{m_bk_benign.group(1)}A"
            if source_id in by_source_id:
                mapping[expr_id] = by_source_id[source_id]
                continue

        # c{digits}A -> M2200{digits zero-padded to 3}A (BioBank cancer -> BC samples)
        m_bk_cancer = re.match(r"^c(\d+)A$", prefix)
        if m_bk_cancer:
            source_id = f"M2200{int(m_bk_cancer.group(1)):03d}A"
            if source_id in by_source_id:
                mapping[expr_id] = by_source_id[source_id]
                continue

        # HR10 special case: HR10_Sample141_R1 (capital S)
        m_hr = re.match(r"^(HR\d+)_[Ss]ample", expr_id)
        if m_hr:
            hr_id = m_hr.group(1)
            if hr_id in by_rna_id:
                mapping[expr_id] = by_rna_id[hr_id]
                continue

        unmatched.append(expr_id)

    return mapping, unmatched


# ---------------------------------------------------------------------------
# 4. Produce clean clinical CSV
# ---------------------------------------------------------------------------

OUTPUT_COLS = [
    "SampleID",           # = expression row name
    "Group",              # Tumor / Normal
    "Age",                # numeric
    "Gender",
    "Histology",
    # Overall stage
    "AJCC_Stage",         # e.g. Stage IIA
    "Stage_Group",        # e.g. Stage II  (for ROC grouping)
    "Stage_Binary",       # Early (0/I/II) vs Late (III/IV)
    # T stage
    "Clinical_T",         # raw standardized
    "Pathologic_T",       # raw standardized
    "T_Stage",            # best available (path > clin)
    "T_Category",         # T0/Tis/T1/T2/T3/T4 (for grouping)
    # N stage
    "Clinical_N",
    "Pathologic_N",
    "N_Stage",
    "N_Category",         # N0/N1/N2/N3
    # M stage
    "Clinical_M",
    "Pathologic_M",
    "M_Stage",
    "M_Category",         # M0/M1
    # Receptor status
    "ER_Status",
    "PR_Status",
    "HER2_Status",
    "Subtype",            # Luminal A / Luminal B / HER2+ / Triple-negative
    # Survival (batch2 clinical only; batch1 has none)
    "OS_time",            # months
    "OS_status",          # 0=censored, 1=event
    "RFS_time",           # recurrence-free survival months
    "Alive",              # original value
    "Recurrence",
    # Treatment
    "Neo_Adjuvant",
    "Adjuvant_Rx",
    # Source info
    "source_ID",
    "RNA_ID",
    "Cohort",
    "Batch",              # batch1 or batch2
]


def build_clean_row(expr_id, clin_row, batch_label):
    """Build a standardized output row from an expression ID + clinical row."""
    out = {col: "" for col in OUTPUT_COLS}
    out["SampleID"] = expr_id

    if clin_row is None:
        # No clinical data for this sample (expression-only)
        out["Group"] = ""
        out["Batch"] = batch_label
        return out

    # Group
    group_raw = clean_str(clin_row.get("group", ""))
    if group_raw.lower() in ("cancer", "tumor"):
        out["Group"] = "Tumor"
    elif group_raw.lower() in ("benign", "normal", "health"):
        out["Group"] = "Normal"
    else:
        out["Group"] = group_raw

    # Age - handle different column names
    age_raw = clean_str(
        clin_row.get("Age at Sample Acquisition", "") or
        clin_row.get("Age at Collection", "")
    )
    try:
        out["Age"] = int(float(age_raw))
    except (ValueError, TypeError):
        out["Age"] = ""

    out["Gender"] = clean_str(
        clin_row.get("Gender", "")
    )
    out["Histology"] = clean_str(
        clin_row.get("Histologic Type ( First letter Capitalize) ", "") or
        clin_row.get("Histology", "")
    )

    # --- Staging ---
    # Overall AJCC stage
    stage_raw = clean_str(clin_row.get("Best CS/AJCC Stage", ""))
    ajcc_stage, stage_group = parse_overall_stage(stage_raw)
    out["AJCC_Stage"] = ajcc_stage
    out["Stage_Group"] = stage_group

    # Stage binary: Early (0/I/II) vs Late (III/IV)
    if stage_group:
        if stage_group in ("Stage 0", "Stage I", "Stage II"):
            out["Stage_Binary"] = "Early"
        elif stage_group in ("Stage III", "Stage IV"):
            out["Stage_Binary"] = "Late"
        else:
            out["Stage_Binary"] = ""

    # Clinical T/N/M
    clin_t = clean_str(clin_row.get("AJCC Clinical T", ""))
    clin_n = clean_str(clin_row.get("AJCC Clinical N", ""))
    clin_m = clean_str(clin_row.get("AJCC Clinical M", ""))

    # Pathologic T/N/M - handle both column naming conventions
    path_t = clean_str(
        clin_row.get("AJCC Pathologic T", "") or
        clin_row.get("T Stage", "")
    )
    path_n = clean_str(
        clin_row.get("AJCC Pathologic N", "") or
        clin_row.get("N Stage", "")
    )
    path_m = clean_str(
        clin_row.get("AJCC Pathologic M", "") or
        clin_row.get("M stage", "")
    )

    # Parse and standardize
    ct_val, _ = parse_tnm(clin_t)
    pt_val, _ = parse_tnm(path_t)
    out["Clinical_T"] = f"cT{ct_val}" if ct_val else ""
    out["Pathologic_T"] = f"pT{pt_val}" if pt_val else ""
    t_stage, t_cat = standardize_t_stage(clin_t, path_t)
    out["T_Stage"] = t_stage
    out["T_Category"] = t_cat

    cn_val, _ = parse_tnm(clin_n)
    pn_val, _ = parse_tnm(path_n)
    out["Clinical_N"] = f"cN{cn_val}" if cn_val else ""
    out["Pathologic_N"] = f"pN{pn_val}" if pn_val else ""
    n_stage, n_cat = standardize_n_stage(clin_n, path_n)
    out["N_Stage"] = n_stage
    out["N_Category"] = n_cat

    cm_val, _ = parse_tnm(clin_m)
    pm_val, _ = parse_tnm(path_m)
    out["Clinical_M"] = f"cM{cm_val}" if cm_val else ""
    out["Pathologic_M"] = f"pM{pm_val}" if pm_val else ""
    m_stage, m_cat = standardize_m_stage(clin_m, path_m)
    out["M_Stage"] = m_stage
    out["M_Category"] = m_cat

    # --- Receptor status & subtype ---
    er_raw = clean_str(
        clin_row.get("ER ", "") or clin_row.get("ER Status", "")
    )
    pr_raw = clean_str(
        clin_row.get("PR  ", "") or clin_row.get("PR Status", "")
    )
    her2_raw = clean_str(
        clin_row.get("Her2                                                 (IHC)", "") or
        clin_row.get("HER2 Status", "")
    )
    subtype, er, pr, her2 = derive_subtype(er_raw, pr_raw, her2_raw)
    out["ER_Status"] = er
    out["PR_Status"] = pr
    out["HER2_Status"] = her2
    out["Subtype"] = subtype if out["Group"] == "Tumor" else ""

    # --- Survival (only batch2 clinical has these) ---
    alive_raw = clean_str(clin_row.get("Alive Yes/No", ""))
    os_mo_raw = clean_str(clin_row.get("OS (Mo)", ""))
    rfs_raw = clean_str(clin_row.get("Recurrence Free Survival (Mo)", ""))
    recur_raw = clean_str(clin_row.get("Reccurence/Progression", ""))

    out["Alive"] = alive_raw
    out["OS_time"] = parse_os_months(os_mo_raw)
    os_status = parse_alive_status(alive_raw)
    out["OS_status"] = os_status
    out["RFS_time"] = parse_os_months(rfs_raw)
    out["Recurrence"] = recur_raw

    # Treatment
    out["Neo_Adjuvant"] = clean_str(clin_row.get("Neo-AdjuvantTreatment", ""))
    out["Adjuvant_Rx"] = clean_str(clin_row.get("Adjuvant Rx", ""))

    # Source info
    out["source_ID"] = clean_str(
        clin_row.get("source_ID", "") or clin_row.get("source ID", "")
    )
    out["RNA_ID"] = clean_str(clin_row.get("RNA_ID", ""))
    out["Cohort"] = clean_str(clin_row.get("Cohort", ""))
    out["Batch"] = batch_label

    return out


def write_clean_csv(rows, path):
    """Write clean CSV."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLS)
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# 5. Main
# ---------------------------------------------------------------------------

def main():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Read clinical data
    batch1_clin = read_csv(os.path.join(base, "clinical_data", "yyfbatch1_clinical.csv"))
    batch2_clin = read_csv(os.path.join(base, "clinical_data", "yyfbatch2_clinical.csv"))

    # Read expression IDs
    batch1_expr_ids = get_expression_ids(
        os.path.join(base, "processed_results", "yyfbatch1_processed.csv"))
    batch2_expr_ids = get_expression_ids(
        os.path.join(base, "processed_results", "yyfbatch2_processed.csv"))

    print(f"Batch1 expression samples: {len(batch1_expr_ids)}")
    print(f"Batch2 expression samples: {len(batch2_expr_ids)}")
    print(f"Batch1 clinical rows: {len(batch1_clin)}")
    print(f"Batch2 clinical rows: {len(batch2_clin)}")
    print()

    # -----------------------------------------------------------------------
    # Batch1 expression <-> Batch2 clinical (via RNA_ID2)
    # -----------------------------------------------------------------------
    b2c_map = build_batch1expr_to_batch2clin(batch2_clin)
    batch1_clean = []
    batch1_matched = 0
    for eid in batch1_expr_ids:
        clin_row = b2c_map.get(eid)
        if clin_row:
            batch1_matched += 1
        batch1_clean.append(build_clean_row(eid, clin_row, "yyfbatch1"))

    print(f"Batch1 (yyfbatch1_processed): {batch1_matched}/{len(batch1_expr_ids)} "
          f"matched to batch2 clinical")

    # -----------------------------------------------------------------------
    # Batch2 expression <-> Batch1 clinical (via RNA_ID prefix)
    # -----------------------------------------------------------------------
    b1c_map, unmatched = build_batch2expr_to_batch1clin(batch1_clin, batch2_expr_ids)
    batch2_clean = []
    batch2_matched = 0
    for eid in batch2_expr_ids:
        clin_row = b1c_map.get(eid)
        if clin_row:
            batch2_matched += 1
        batch2_clean.append(build_clean_row(eid, clin_row, "yyfbatch2"))

    print(f"Batch2 (yyfbatch2_processed): {batch2_matched}/{len(batch2_expr_ids)} "
          f"matched to batch1 clinical")
    if unmatched:
        print(f"  Unmatched batch2 expression IDs ({len(unmatched)}): {unmatched[:10]}...")

    # -----------------------------------------------------------------------
    # Write output
    # -----------------------------------------------------------------------
    out_dir = os.path.join(base, "clinical_data")
    write_clean_csv(batch1_clean,
                    os.path.join(out_dir, "yyfbatch1_clinical_clean.csv"))
    write_clean_csv(batch2_clean,
                    os.path.join(out_dir, "yyfbatch2_clinical_clean.csv"))

    print(f"\nOutput written:")
    print(f"  clinical_data/yyfbatch1_clinical_clean.csv ({len(batch1_clean)} rows)")
    print(f"  clinical_data/yyfbatch2_clinical_clean.csv ({len(batch2_clean)} rows)")

    # -----------------------------------------------------------------------
    # Summary stats
    # -----------------------------------------------------------------------
    for label, rows in [("yyfbatch1", batch1_clean), ("yyfbatch2", batch2_clean)]:
        tumors = [r for r in rows if r["Group"] == "Tumor"]
        normals = [r for r in rows if r["Group"] == "Normal"]
        staged = [r for r in tumors if r["AJCC_Stage"]]
        subtypes = {}
        for r in tumors:
            st = r["Subtype"]
            if st:
                subtypes[st] = subtypes.get(st, 0) + 1
        stage_groups = {}
        for r in tumors:
            sg = r["Stage_Group"]
            if sg:
                stage_groups[sg] = stage_groups.get(sg, 0) + 1
        t_cats = {}
        for r in tumors:
            tc = r["T_Category"]
            if tc:
                t_cats[tc] = t_cats.get(tc, 0) + 1
        n_cats = {}
        for r in tumors:
            nc = r["N_Category"]
            if nc:
                n_cats[nc] = n_cats.get(nc, 0) + 1
        has_os = sum(1 for r in tumors if r["OS_time"] != "")

        print(f"\n{'='*60}")
        print(f"  {label}: {len(tumors)} Tumor, {len(normals)} Normal")
        print(f"  Staged: {len(staged)}/{len(tumors)} tumors")
        print(f"  Stage groups: {stage_groups}")
        print(f"  T categories: {t_cats}")
        print(f"  N categories: {n_cats}")
        print(f"  Subtypes: {subtypes}")
        print(f"  Has OS data: {has_os}/{len(tumors)} tumors")


if __name__ == "__main__":
    main()
