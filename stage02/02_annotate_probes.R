# ============================================================
# Stage 02: Probe Annotation
# Dataset: GSE45827
# Platform: Affymetrix Human Genome U133 Plus 2.0
# Annotation: hgu133plus2.db
# ============================================================

source("renv/activate.R")

library(hgu133plus2.db)
library(AnnotationDbi)

results_dir <- "stage02/results"

contrast_files <- list.files(
  results_dir,
  pattern = "^DE_(Basal_vs_Her2|Basal_vs_LuminalA|Basal_vs_LuminalB|Her2_vs_LuminalA|Her2_vs_LuminalB|LuminalB_vs_LuminalA)\\.csv$",
  full.names = TRUE
)

# Exclude summary table
contrast_files <- contrast_files[
  !grepl("DE_summary\\.csv$", contrast_files)
]

cat("Files to annotate:\n")
print(contrast_files)

# ------------------------------------------------------------
# Annotation function
# ------------------------------------------------------------

annotate_de <- function(file) {

  cat("\n========================================\n")
  cat("Annotating:", file, "\n")
  cat("========================================\n")

  results <- read.csv(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  stopifnot("Probe_ID" %in% colnames(results))

  annotation <- AnnotationDbi::select(
    hgu133plus2.db,
    keys = unique(results$Probe_ID),
    columns = c(
      "SYMBOL",
      "ENTREZID",
      "ENSEMBL",
      "GENENAME"
    ),
    keytype = "PROBEID"
  )

  cat("DE rows:", nrow(results), "\n")
  cat("Annotation rows:", nrow(annotation), "\n")

  # ----------------------------------------------------------
  # Mapping diagnostics
  # ----------------------------------------------------------

  input_probes <- unique(results$Probe_ID)
  mapped_probes <- unique(annotation$PROBEID)

  n_mapped <- sum(input_probes %in% mapped_probes)
  n_unmapped <- sum(!input_probes %in% mapped_probes)

  cat("Mapped probes:", n_mapped, "\n")
  cat("Unmapped probes:", n_unmapped, "\n")

  # Number of probes with multiple annotation rows
  annotation_counts <- table(annotation$PROBEID)

  multi_mapping <- names(
    annotation_counts[annotation_counts > 1]
  )

  n_multi <- length(multi_mapping)

  n_unique <- n_mapped - n_multi

  cat(
    "Uniquely annotated probes:",
    n_unique,
    "\n"
  )

  cat(
    "Probes with multiple annotation rows:",
    n_multi,
    "\n"
  )

  # ----------------------------------------------------------
  # Collapse multiple annotation rows
  # ----------------------------------------------------------

  annotation_collapsed <- aggregate(
    cbind(
      SYMBOL,
      ENTREZID,
      ENSEMBL,
      GENENAME
    ) ~ PROBEID,
    data = annotation,
    FUN = function(x) {

      x <- unique(
        x[!is.na(x) & x != ""]
      )

      if (length(x) == 0) {
        return(NA_character_)
      }

      paste(x, collapse = ";")
    }
  )

  # ----------------------------------------------------------
  # Merge annotation with DE results
  # ----------------------------------------------------------

  annotated <- merge(
    results,
    annotation_collapsed,
    by.x = "Probe_ID",
    by.y = "PROBEID",
    all.x = TRUE,
    sort = FALSE
  )

  # Restore original DE ordering
  annotated <- annotated[
    match(results$Probe_ID, annotated$Probe_ID),
  ]

  rownames(annotated) <- NULL

  # ----------------------------------------------------------
  # Save
  # ----------------------------------------------------------

  output_file <- sub(
    "\\.csv$",
    "_annotated.csv",
    file
  )

  write.csv(
    annotated,
    output_file,
    row.names = FALSE
  )

  cat("Saved:", output_file, "\n")

  return(annotated)
}

# ------------------------------------------------------------
# Run annotation for all contrasts
# ------------------------------------------------------------

annotated_results <- lapply(
  contrast_files,
  annotate_de
)

names(annotated_results) <- basename(
  contrast_files
)

cat("\n========================================\n")
cat("Probe annotation complete\n")
cat("========================================\n")