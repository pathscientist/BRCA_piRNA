################################################################################
#                                                                              #
#   Reformat TCGA-BRCA Clinical Data                                          #
#                                                                              #
#   Reads:                                                                    #
#     1. clinical_data/TCGA_BRCA_clinical.csv  (raw TCGA column names)        #
#     2. clinical_data/clinical.tsv            (GDC portal download for Age)  #
#                                                                              #
#   Produces:                                                                 #
#     clinical_data/TCGA_BRCA_clinical.csv  (pipeline-ready format)           #
#     Columns: SampleID, PatientID, Group, Batch, Age, Gender, Stage,         #
#              Stage_Raw, Stage_Binary, ER_Status, PR_Status, HER2_Status,    #
#              Subtype, OS_time, OS_status                                    #
#                                                                              #
#   Usage: Rscript scripts/reformat_tcga_clinical.R                           #
#                                                                              #
################################################################################

cat("=== Reformat TCGA-BRCA Clinical Data ===\n\n")

# ==============================================================================
# 1. READ RAW CLINICAL CSV
# ==============================================================================

raw_path <- "clinical_data/TCGA_BRCA_clinical.csv"
if (!file.exists(raw_path)) stop("Raw clinical file not found: ", raw_path)

# Read raw CSV — has unnamed/empty columns and duplicate column names
raw <- read.csv(raw_path, stringsAsFactors = FALSE, check.names = FALSE,
                na.strings = c("", "NA", "[Not Available]", "[Not Applicable]",
                               "not reported", "Not Reported"))

cat(sprintf("Loaded raw clinical CSV: %d rows, %d columns\n", nrow(raw), ncol(raw)))

# Assign usable column names (raw has empty names and duplicates)
colnames(raw) <- c(
  "SampleID", "SampleID2",
  "ER_Status_raw", "PR_Status_raw",
  "HER2_IHC_level", "HER2_IHC_status",
  "HER2_FISH", "HER2_merge",
  "X1",  # unnamed column

  "tumor_stage_raw",
  "days_to_new_tumor_event",
  "days_to_death",
  "days_to_last_follow_up",
  "X2", "X3",  # unnamed columns
  "person_neoplasm_cancer_status",
  "vital_status",
  "X4"  # trailing column
)

cat(sprintf("  Sample IDs: %s ... %s\n", raw$SampleID[1], raw$SampleID[nrow(raw)]))

# ==============================================================================
# 2. READ GDC clinical.tsv (for Age, gender, etc.)
# ==============================================================================

gdc_path <- "clinical_data/clinical.tsv"
gdc <- NULL

if (file.exists(gdc_path)) {
  gdc <- read.delim(gdc_path, stringsAsFactors = FALSE,
                    na.strings = c("", "NA", "'--", "not reported", "Not Reported"))
  cat(sprintf("Loaded GDC clinical.tsv: %d rows, %d columns\n", nrow(gdc), ncol(gdc)))

  # Show available columns for debugging
  age_candidates <- grep("age", colnames(gdc), ignore.case = TRUE, value = TRUE)
  cat(sprintf("  Age columns found: %s\n", paste(age_candidates, collapse = ", ")))
} else {
  cat("WARNING: clinical_data/clinical.tsv not found.\n")
  cat("  Age will be NA. To fix, download from GDC Data Portal:\n")
  cat("  https://portal.gdc.cancer.gov/projects/TCGA-BRCA → Clinical tab → TSV\n")
  cat("  Save as clinical_data/clinical.tsv\n\n")
}

# ==============================================================================
# 3. BUILD STANDARDIZED DATA FRAME
# ==============================================================================

cat("\nBuilding standardized clinical table...\n")

clin <- data.frame(
  SampleID  = raw$SampleID,
  PatientID = substr(raw$SampleID, 1, 12),
  stringsAsFactors = FALSE
)

# --- Group (Tumor vs Normal) from TCGA barcode ---
sample_type_code <- as.numeric(substr(clin$SampleID, 14, 15))
clin$Group <- ifelse(sample_type_code >= 1 & sample_type_code <= 9, "Tumor",
              ifelse(sample_type_code >= 10 & sample_type_code <= 19, "Normal", NA))
clin$Batch <- "BRCA1"

# --- Age from GDC clinical.tsv ---
clin$Age <- NA_real_
if (!is.null(gdc)) {
  # Try common GDC column names for age
  age_col <- NULL
  for (ac in c("age_at_diagnosis", "age_at_index",
               "age_at_initial_pathologic_diagnosis")) {
    if (ac %in% colnames(gdc)) { age_col <- ac; break }
  }

  if (!is.null(age_col)) {
    age_vals <- as.numeric(gdc[[age_col]])

    # GDC age_at_diagnosis is in days → convert to years
    if (median(age_vals, na.rm = TRUE) > 200) {
      age_vals <- round(age_vals / 365.25)
      cat(sprintf("  Age source: %s (converted from days to years)\n", age_col))
    } else {
      cat(sprintf("  Age source: %s (already in years)\n", age_col))
    }

    # Build patient → age lookup
    id_col <- NULL
    for (ic in c("submitter_id", "case_submitter_id", "bcr_patient_barcode")) {
      if (ic %in% colnames(gdc)) { id_col <- ic; break }
    }

    if (!is.null(id_col)) {
      age_lookup <- setNames(age_vals, gdc[[id_col]])
      # Remove duplicates (keep first)
      age_lookup <- age_lookup[!duplicated(names(age_lookup))]
      matched_ages <- age_lookup[clin$PatientID]
      clin$Age <- as.numeric(matched_ages)
      cat(sprintf("  Age matched: %d / %d patients\n",
                  sum(!is.na(clin$Age)), nrow(clin)))
    }
  }

  # --- Gender from GDC ---
  clin$Gender <- NA_character_
  for (gc in c("gender", "demographic.gender")) {
    if (gc %in% colnames(gdc)) {
      gender_lookup <- setNames(tolower(gdc[[gc]]), gdc[[id_col]])
      gender_lookup <- gender_lookup[!duplicated(names(gender_lookup))]
      clin$Gender <- as.character(gender_lookup[clin$PatientID])
      break
    }
  }
}

# --- ER status ---
clin$ER_Status <- tolower(trimws(raw$ER_Status_raw))
clin$ER_Status <- ifelse(clin$ER_Status %in% c("positive", "pos"), "positive",
                  ifelse(clin$ER_Status %in% c("negative", "neg"), "negative", NA))

# --- PR status ---
clin$PR_Status <- tolower(trimws(raw$PR_Status_raw))
clin$PR_Status <- ifelse(clin$PR_Status %in% c("positive", "pos"), "positive",
                  ifelse(clin$PR_Status %in% c("negative", "neg"), "negative", NA))

# --- HER2 status (use the merged call, fall back to IHC status) ---
her2_raw <- tolower(trimws(raw$HER2_merge))
her2_raw <- ifelse(is.na(her2_raw), tolower(trimws(raw$HER2_IHC_status)), her2_raw)
clin$HER2_Status <- ifelse(her2_raw %in% c("positive", "pos"), "positive",
                   ifelse(her2_raw %in% c("negative", "neg"), "negative",
                   ifelse(her2_raw %in% c("equivocal", "indeterminate"), "equivocal", NA)))

# --- Stage ---
stage_raw <- tolower(trimws(raw$tumor_stage_raw))
clin$Stage_Raw <- raw$tumor_stage_raw

clin$Stage <- sapply(stage_raw, function(s) {
  if (is.na(s)) return(NA)
  s <- gsub("^stage\\s*", "", s)
  if (grepl("^iv", s))   return("Stage IV")
  if (grepl("^iii", s))  return("Stage III")
  if (grepl("^ii", s))   return("Stage II")
  if (grepl("^i", s))    return("Stage I")
  if (grepl("^x", s))    return(NA)  # stage x = not available
  return(NA)
})

clin$Stage_Binary <- ifelse(clin$Stage %in% c("Stage III", "Stage IV"), "Late", "Early")

# --- Molecular Subtype ---
clin$Subtype <- mapply(function(er, pr, her2) {
  if (is.na(er) && is.na(pr) && is.na(her2)) return(NA)
  er_pos  <- !is.na(er)  && er  == "positive"
  pr_pos  <- !is.na(pr)  && pr  == "positive"
  her2_pos <- !is.na(her2) && her2 == "positive"
  her2_neg <- !is.na(her2) && her2 == "negative"
  all_neg <- !er_pos && !pr_pos && her2_neg

  if (all_neg) return("Triple-negative")
  if (her2_pos && (er_pos || pr_pos)) return("Luminal B")
  if (her2_pos) return("HER2+")
  if (er_pos || pr_pos) return("Luminal A")
  return("Unclassified")
}, clin$ER_Status, clin$PR_Status, clin$HER2_Status, USE.NAMES = FALSE)

# --- Survival: OS_time and OS_status ---
vital <- tolower(trimws(raw$vital_status))
clin$OS_status <- ifelse(vital == "dead", 1,
                  ifelse(vital %in% c("alive", "living"), 0, NA))

dtd  <- as.numeric(raw$days_to_death)
dtlfu <- as.numeric(raw$days_to_last_follow_up)

# For dead → days_to_death; for alive → days_to_last_follow_up
os_days <- ifelse(clin$OS_status == 1 & !is.na(dtd), dtd, dtlfu)
clin$OS_time <- round(os_days / 30.44, 2)  # convert to months

# ==============================================================================
# 4. SAVE
# ==============================================================================

# Select output columns
output_cols <- c("SampleID", "PatientID", "Group", "Batch",
                 "Age", "Gender", "Stage", "Stage_Raw", "Stage_Binary",
                 "ER_Status", "PR_Status", "HER2_Status", "Subtype",
                 "OS_time", "OS_status")
output_cols <- intersect(output_cols, colnames(clin))

clin_final <- clin[, output_cols]

# Remove duplicates
clin_final <- clin_final[!duplicated(clin_final$SampleID), ]
clin_final <- clin_final[order(clin_final$SampleID), ]

# Write output — use a separate path if the original is read-only
out_path <- raw_path
if (file.exists(raw_path) && file.access(raw_path, 2) != 0) {
  out_path <- sub("\\.csv$", "_reformatted.csv", raw_path)
  cat(sprintf("  Original is read-only — writing to: %s\n", out_path))
}

write.csv(clin_final, out_path, row.names = FALSE)

# ==============================================================================
# 5. SUMMARY
# ==============================================================================

cat(sprintf("\nSaved reformatted clinical data: %s\n", out_path))
cat(sprintf("  Total records: %d\n", nrow(clin_final)))
cat(sprintf("  Tumor:  %d\n", sum(clin_final$Group == "Tumor", na.rm = TRUE)))
cat(sprintf("  Normal: %d\n", sum(clin_final$Group == "Normal", na.rm = TRUE)))

cat(sprintf("\n  Age:     %d patients with data",
            sum(!is.na(clin_final$Age))))
if (any(!is.na(clin_final$Age))) {
  cat(sprintf(" (median=%.0f, range=%d-%d)",
              median(clin_final$Age, na.rm = TRUE),
              min(clin_final$Age, na.rm = TRUE),
              max(clin_final$Age, na.rm = TRUE)))
}

cat("\n\n  Stage:\n")
print(table(clin_final$Stage, useNA = "ifany"))

cat("\n  Subtype:\n")
print(table(clin_final$Subtype, useNA = "ifany"))

cat(sprintf("\n  OS_time:   %d patients with data",
            sum(!is.na(clin_final$OS_time))))
if (any(!is.na(clin_final$OS_time))) {
  cat(sprintf(" (median: %.1f months)", median(clin_final$OS_time, na.rm = TRUE)))
}
cat(sprintf("\n  OS_status: %d events / %d total\n",
            sum(clin_final$OS_status == 1, na.rm = TRUE),
            sum(!is.na(clin_final$OS_status))))

cat("\n=== REFORMATTING COMPLETE ===\n")
cat("\nThe pipeline will now read standardized columns:\n")
cat("  SampleID, Age, Stage, Subtype, OS_time, OS_status\n")
cat("\nTo add Age data, download clinical.tsv from GDC:\n")
cat("  https://portal.gdc.cancer.gov/projects/TCGA-BRCA\n")
cat("  Repository → Clinical tab → TSV download → save as clinical_data/clinical.tsv\n")
cat("  Then re-run: Rscript scripts/reformat_tcga_clinical.R\n")
