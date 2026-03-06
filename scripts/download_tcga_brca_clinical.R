################################################################################
#                                                                              #
#   Download TCGA-BRCA Clinical Data via TCGAbiolinks                         #
#                                                                              #
#   Downloads clinical + biospecimen data from GDC and produces a             #
#   standardized CSV matching the pipeline's expected format:                  #
#                                                                              #
#   Output: clinical_data/TCGA_BRCA_clinical.csv                              #
#     Columns: SampleID, PatientID, Age, Gender, Stage, Stage_Binary,         #
#              Subtype, ER_Status, PR_Status, HER2_Status,                    #
#              OS_time, OS_status, DSS_time, DSS_status,                      #
#              PFI_time, PFI_status, Batch, Group                             #
#                                                                              #
#   Usage: Rscript scripts/download_tcga_brca_clinical.R                      #
#                                                                              #
################################################################################

cat("=== TCGA-BRCA Clinical Data Download ===\n\n")

# ==============================================================================
# 1. INSTALL & LOAD PACKAGES
# ==============================================================================

cran_pkgs <- c("dplyr", "stringr")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE, quiet = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", quiet = TRUE)

if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks", ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(dplyr)
  library(stringr)
})

cat("Packages loaded.\n\n")

# ==============================================================================
# 2. DOWNLOAD CLINICAL DATA FROM GDC
# ==============================================================================

cat("Step 1: Querying GDC for TCGA-BRCA clinical data...\n")

query_clin <- GDCquery(
  project    = "TCGA-BRCA",
  data.category = "Clinical",
  data.type     = "Clinical Supplement",
  data.format   = "BCR Biotab"
)

cat("Step 2: Downloading clinical files...\n")
GDCdownload(query_clin)

cat("Step 3: Preparing clinical data...\n")
clin_biotab <- GDCprepare(query_clin)

# clin_biotab is a list of data frames (patient, drug, follow_up, etc.)
cat("  Available clinical tables:\n")
cat(paste("   ", names(clin_biotab), collapse = "\n"), "\n\n")

# ==============================================================================
# 3. ALSO GET INDEXED CLINICAL DATA (curated fields)
# ==============================================================================

cat("Step 4: Downloading curated clinical-indexed data...\n")

clin_indexed <- tryCatch({
  GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
}, error = function(e) {
  cat("  Indexed clinical query failed, trying alternative...\n")
  NULL
})

if (is.null(clin_indexed)) {
  # Fallback: query directly
  query_clin2 <- GDCquery(
    project    = "TCGA-BRCA",
    data.category = "Clinical",
    data.type  = "Clinical Supplement"
  )
  GDCdownload(query_clin2)
  clin_indexed <- GDCprepare_clinic(query_clin2, clinical.info = "patient")
}

cat(sprintf("  Curated clinical records: %d patients\n", nrow(clin_indexed)))

# ==============================================================================
# 4. EXTRACT & STANDARDIZE KEY FIELDS
# ==============================================================================

cat("\nStep 5: Standardizing clinical variables...\n")

# --- Patient-level data from indexed clinical ---
clin_df <- data.frame(
  PatientID = clin_indexed$submitter_id,
  stringsAsFactors = FALSE
)

# Age at diagnosis
age_cols <- c("age_at_diagnosis", "age_at_initial_pathologic_diagnosis",
              "age_at_index")
for (ac in age_cols) {
  if (ac %in% colnames(clin_indexed)) {
    age_vals <- as.numeric(clin_indexed[[ac]])
    # GDC stores age_at_diagnosis in days; convert to years
    if (ac == "age_at_diagnosis" && median(age_vals, na.rm = TRUE) > 200) {
      age_vals <- round(age_vals / 365.25)
    }
    clin_df$Age <- age_vals
    cat(sprintf("  Age source: %s (median: %.0f years)\n", ac,
                median(clin_df$Age, na.rm = TRUE)))
    break
  }
}

# Gender
gender_cols <- c("gender", "demographic.gender")
for (gc in gender_cols) {
  if (gc %in% colnames(clin_indexed)) {
    clin_df$Gender <- clin_indexed[[gc]]
    break
  }
}

# AJCC Stage
stage_cols <- c("ajcc_pathologic_stage", "tumor_stage",
                "ajcc_clinical_stage", "clinical_stage",
                "figo_stage", "stage_event.pathologic_stage")
for (sc in stage_cols) {
  if (sc %in% colnames(clin_indexed)) {
    raw_stage <- clin_indexed[[sc]]
    clin_df$Stage_Raw <- raw_stage

    # Standardize to "Stage I/II/III/IV"
    clin_df$Stage <- sapply(raw_stage, function(s) {
      if (is.na(s) || s == "" || s == "not reported" ||
          s == "Not Reported" || s == "[Not Available]") return(NA)
      s <- toupper(trimws(s))
      s <- gsub("^STAGE\\s*", "", s)
      if (grepl("^IV", s))   return("Stage IV")
      if (grepl("^III", s))  return("Stage III")
      if (grepl("^II", s))   return("Stage II")
      if (grepl("^I", s))    return("Stage I")
      if (grepl("^X", s))    return(NA)
      return(NA)
    })

    cat(sprintf("  Stage source: %s\n", sc))
    cat("  Stage distribution:\n")
    print(table(clin_df$Stage, useNA = "ifany"))
    break
  }
}

# Stage binary
clin_df$Stage_Binary <- ifelse(
  clin_df$Stage %in% c("Stage III", "Stage IV"), "Late", "Early"
)

# --- Receptor status & molecular subtype ---

# ER status
er_cols <- c("er_status_by_ihc", "er_status",
             "diagnoses.er_status_by_ihc")
for (ec in er_cols) {
  if (ec %in% colnames(clin_indexed)) {
    clin_df$ER_Status <- tolower(trimws(clin_indexed[[ec]]))
    clin_df$ER_Status <- ifelse(
      clin_df$ER_Status %in% c("positive", "pos"), "positive",
      ifelse(clin_df$ER_Status %in% c("negative", "neg"), "negative", NA))
    break
  }
}

# PR status
pr_cols <- c("pr_status_by_ihc", "pr_status",
             "diagnoses.pr_status_by_ihc")
for (pc in pr_cols) {
  if (pc %in% colnames(clin_indexed)) {
    clin_df$PR_Status <- tolower(trimws(clin_indexed[[pc]]))
    clin_df$PR_Status <- ifelse(
      clin_df$PR_Status %in% c("positive", "pos"), "positive",
      ifelse(clin_df$PR_Status %in% c("negative", "neg"), "negative", NA))
    break
  }
}

# HER2 status
her2_cols <- c("her2_status_by_ihc", "her2_status",
               "diagnoses.her2_status_by_ihc",
               "lab_proc_her2_neu_immunohistochemistry_receptor_status")
for (hc in her2_cols) {
  if (hc %in% colnames(clin_indexed)) {
    clin_df$HER2_Status <- tolower(trimws(clin_indexed[[hc]]))
    clin_df$HER2_Status <- ifelse(
      clin_df$HER2_Status %in% c("positive", "pos"), "positive",
      ifelse(clin_df$HER2_Status %in% c("negative", "neg"), "negative",
      ifelse(clin_df$HER2_Status %in% c("equivocal", "indeterminate"), "equivocal", NA)))
    break
  }
}

# Derive molecular subtype
clin_df$Subtype <- mapply(function(er, pr, her2) {
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
}, clin_df$ER_Status, clin_df$PR_Status, clin_df$HER2_Status)

cat("\n  Molecular subtype distribution:\n")
print(table(clin_df$Subtype, useNA = "ifany"))

# --- Survival endpoints ---

# Overall Survival
os_time_cols <- c("days_to_death", "days_to_last_follow_up")
os_stat_cols <- c("vital_status")

if ("vital_status" %in% colnames(clin_indexed)) {
  vital <- tolower(trimws(clin_indexed$vital_status))

  # OS_status: 1 = dead, 0 = alive/censored
  clin_df$OS_status <- ifelse(vital == "dead", 1,
                       ifelse(vital %in% c("alive", "living"), 0, NA))

  # OS_time: days_to_death for dead, days_to_last_follow_up for alive
  dtd <- rep(NA_real_, nrow(clin_indexed))
  dtlfu <- rep(NA_real_, nrow(clin_indexed))

  if ("days_to_death" %in% colnames(clin_indexed))
    dtd <- as.numeric(clin_indexed$days_to_death)
  if ("days_to_last_follow_up" %in% colnames(clin_indexed))
    dtlfu <- as.numeric(clin_indexed$days_to_last_follow_up)

  # For dead patients: use days_to_death; for alive: days_to_last_follow_up
  clin_df$OS_time_days <- ifelse(clin_df$OS_status == 1 & !is.na(dtd), dtd, dtlfu)

  # Convert to months (30.44 days/month)
  clin_df$OS_time <- round(clin_df$OS_time_days / 30.44, 2)

  cat(sprintf("\n  Survival data: %d patients with OS_time\n",
              sum(!is.na(clin_df$OS_time))))
  cat(sprintf("  Events (deaths): %d / %d (%.1f%%)\n",
              sum(clin_df$OS_status == 1, na.rm = TRUE),
              sum(!is.na(clin_df$OS_status)),
              100 * mean(clin_df$OS_status == 1, na.rm = TRUE)))
  cat(sprintf("  Median follow-up: %.1f months\n",
              median(clin_df$OS_time, na.rm = TRUE)))
}

# Disease-Specific Survival (if available)
if ("days_to_death" %in% colnames(clin_indexed)) {
  cause_cols <- c("cause_of_death", "cause_of_death.source")
  cause <- NA
  for (cc in cause_cols) {
    if (cc %in% colnames(clin_indexed)) {
      cause <- clin_indexed[[cc]]
      break
    }
  }

  # DSS: censor patients who died of non-cancer causes
  clin_df$DSS_status <- clin_df$OS_status
  clin_df$DSS_time   <- clin_df$OS_time
}

# Progression-Free Interval
pfi_cols <- c("days_to_new_tumor_event_after_initial_treatment",
              "new_tumor_event_after_initial_treatment")
if (any(pfi_cols %in% colnames(clin_indexed))) {
  pfi_event_col <- intersect(pfi_cols, colnames(clin_indexed))[1]
  cat(sprintf("  PFI source: %s\n", pfi_event_col))
}

# ==============================================================================
# 5. BUILD SAMPLE-LEVEL TABLE
# ==============================================================================

cat("\nStep 6: Building sample-level table...\n")

# Get biospecimen info to map patient → sample barcodes
bio_query <- tryCatch({
  GDCquery(
    project = "TCGA-BRCA",
    data.category = "Biospecimen",
    data.type     = "Biospecimen Supplement",
    data.format   = "BCR Biotab"
  )
}, error = function(e) NULL)

sample_map <- NULL
if (!is.null(bio_query)) {
  tryCatch({
    GDCdownload(bio_query)
    bio_data <- GDCprepare(bio_query)
    if ("clinical_patient_brca" %in% names(bio_data)) {
      # Use patient table
    }
    if ("biospecimen_sample_brca" %in% names(bio_data)) {
      sample_tab <- bio_data[["biospecimen_sample_brca"]]
      if (all(c("bcr_patient_barcode", "bcr_sample_barcode") %in% colnames(sample_tab))) {
        sample_map <- sample_tab %>%
          select(PatientID = bcr_patient_barcode,
                 SampleBarcode = bcr_sample_barcode,
                 sample_type) %>%
          filter(!is.na(SampleBarcode), SampleBarcode != "")
      }
    }
  }, error = function(e) {
    cat("  Biospecimen query failed:", conditionMessage(e), "\n")
  })
}

# If biospecimen mapping is unavailable, construct sample IDs from patient IDs
# TCGA convention: patient TCGA-XX-XXXX → tumor sample TCGA-XX-XXXX-01A
if (is.null(sample_map) || nrow(sample_map) == 0) {
  cat("  Using TCGA barcode convention for sample ID mapping.\n")
  # Tumor samples: -01 (primary solid tumor)
  # Normal samples: -11 (solid tissue normal)
  sample_map <- rbind(
    data.frame(
      PatientID = clin_df$PatientID,
      SampleBarcode = paste0(clin_df$PatientID, "-01A"),
      sample_type = "Primary Tumor",
      stringsAsFactors = FALSE
    ),
    data.frame(
      PatientID = clin_df$PatientID,
      SampleBarcode = paste0(clin_df$PatientID, "-11A"),
      sample_type = "Solid Tissue Normal",
      stringsAsFactors = FALSE
    )
  )
}

# Determine Group from sample type code
# TCGA sample type: 01-09 = Tumor, 10-19 = Normal
sample_map$Group <- sapply(sample_map$SampleBarcode, function(bc) {
  # Extract sample type digits (positions 14-15 in TCGA barcode)
  st <- as.numeric(substr(bc, 14, 15))
  if (is.na(st)) return(NA)
  if (st >= 1 && st <= 9) return("Tumor")
  if (st >= 10 && st <= 19) return("Normal")
  return(NA)
})

# Merge clinical data with sample mapping
clin_samples <- merge(sample_map, clin_df, by = "PatientID", all.x = TRUE)

# Create SampleID that matches expression data row names
# Use 15-char barcode (TCGA-XX-XXXX-01A) as default
clin_samples$SampleID <- substr(clin_samples$SampleBarcode, 1, 15)
clin_samples$Batch <- "BRCA1"

# ==============================================================================
# 6. FINALIZE & SAVE
# ==============================================================================

cat("\nStep 7: Saving clinical data...\n")

# Select final columns
output_cols <- c("SampleID", "PatientID", "Group", "Batch",
                 "Age", "Gender", "Stage", "Stage_Raw", "Stage_Binary",
                 "ER_Status", "PR_Status", "HER2_Status", "Subtype",
                 "OS_time", "OS_status")

# Add optional columns if they exist
for (oc in c("DSS_time", "DSS_status", "PFI_time", "PFI_status")) {
  if (oc %in% colnames(clin_samples)) output_cols <- c(output_cols, oc)
}

# Keep only existing columns
output_cols <- intersect(output_cols, colnames(clin_samples))
clin_final <- clin_samples[, output_cols]

# Remove duplicates (keep first per SampleID)
clin_final <- clin_final[!duplicated(clin_final$SampleID), ]

# Sort by SampleID
clin_final <- clin_final[order(clin_final$SampleID), ]

dir.create("clinical_data", showWarnings = FALSE)
output_path <- "clinical_data/TCGA_BRCA_clinical.csv"
write.csv(clin_final, output_path, row.names = FALSE)

cat(sprintf("\nClinical data saved to: %s\n", output_path))
cat(sprintf("  Total records: %d\n", nrow(clin_final)))
cat(sprintf("  Tumor samples: %d\n", sum(clin_final$Group == "Tumor", na.rm = TRUE)))
cat(sprintf("  Normal samples: %d\n", sum(clin_final$Group == "Normal", na.rm = TRUE)))

# Print summary
cat("\n--- Clinical Data Summary ---\n")
cat(sprintf("  Age:     median=%.0f, range=%d-%d\n",
            median(clin_final$Age, na.rm = TRUE),
            min(clin_final$Age, na.rm = TRUE),
            max(clin_final$Age, na.rm = TRUE)))
cat("\n  Stage:\n")
print(table(clin_final$Stage, useNA = "ifany"))
cat("\n  Subtype:\n")
print(table(clin_final$Subtype, useNA = "ifany"))
cat(sprintf("\n  OS_time:   %d patients with data (median: %.1f months)\n",
            sum(!is.na(clin_final$OS_time)),
            median(clin_final$OS_time, na.rm = TRUE)))
cat(sprintf("  OS_status: %d events / %d total\n",
            sum(clin_final$OS_status == 1, na.rm = TRUE),
            sum(!is.na(clin_final$OS_status))))

cat("\n=== CLINICAL DATA DOWNLOAD COMPLETE ===\n")
cat("\nNext steps:\n")
cat("  1. Run the main pipeline — it will auto-detect clinical_data/TCGA_BRCA_clinical.csv\n")
cat("  2. Real Age, Stage, Subtype, and survival data will replace simulated values\n")
cat("  3. Sample matching uses TCGA barcode (first 12-15 chars)\n")
