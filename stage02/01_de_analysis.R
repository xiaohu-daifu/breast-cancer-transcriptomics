# ============================================================
# Stage 02: Differential Expression Analysis
# Dataset: GSE45827
# Method: limma
# ============================================================

# -----------------------------
# 1. Load packages
# -----------------------------

library(limma)

# -----------------------------
# 2. File paths
# -----------------------------

expression_file <- "data/GSE45827_expression.csv"
metadata_file <- "data/GSE45827_metadata.csv"

results_dir <- "stage02/results"

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# -----------------------------
# 3. Load data
# -----------------------------

expression <- read.csv(
  expression_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Expression dimensions:", dim(expression), "\n")
cat("Metadata dimensions:", dim(metadata), "\n")

# -----------------------------
# 4. Prepare expression matrix
# -----------------------------

# ID_REF contains Affymetrix probe IDs
probe_ids <- expression$ID_REF

expr <- expression[, -1]

# Convert to numeric matrix
expr <- as.matrix(expr)
mode(expr) <- "numeric"

rownames(expr) <- probe_ids

cat("Expression matrix:", dim(expr), "\n")

# -----------------------------
# 5. Transpose to samples x probes
# -----------------------------

expr <- t(expr)

cat("Transposed expression matrix:", dim(expr), "\n")

# -----------------------------
# 6. Align metadata
# -----------------------------

rownames(metadata) <- metadata$sample_id

metadata <- metadata[rownames(expr), , drop = FALSE]

# Verify sample alignment
stopifnot(identical(
  rownames(expr),
  rownames(metadata)
))

cat("Sample alignment verified.\n")

# -----------------------------
# 7. Define subtype factor
# -----------------------------

metadata$subtype <- factor(
  metadata$subtype,
  levels = c(
    "Luminal A",
    "Luminal B",
    "Her2",
    "Basal"
  )
)

cat("\nSubtype counts:\n")
print(table(metadata$subtype))

# -----------------------------
# 8. Remove probes with low variance
# -----------------------------

probe_variance <- apply(expr, 2, var)

keep <- probe_variance > 0.1

expr_filtered <- expr[, keep, drop = FALSE]

cat("\nProbes before filtering:", ncol(expr), "\n")
cat("Probes after filtering:", ncol(expr_filtered), "\n")
cat("Probes removed:", sum(!keep), "\n")

# -----------------------------
# 9. Create design matrix
# -----------------------------

design <- model.matrix(
  ~ 0 + subtype,
  data = metadata
)

colnames(design) <- c(
  "Luminal_A",
  "Luminal_B",
  "Her2",
  "Basal"
)

cat("\nDesign matrix:\n")
print(head(design))

# -----------------------------
# 10. Fit limma model
# -----------------------------

fit <- lmFit(
  t(expr_filtered),
  design
)

# -----------------------------
# 11. Define biological contrasts
# -----------------------------

contrasts <- makeContrasts(
  Basal_vs_LuminalA = Basal - Luminal_A,
  LuminalB_vs_LuminalA = Luminal_B - Luminal_A,
  Her2_vs_LuminalA = Her2 - Luminal_A,
  Basal_vs_LuminalB = Basal - Luminal_B,
  Her2_vs_LuminalB = Her2 - Luminal_B,
  Basal_vs_Her2 = Basal - Her2,
  levels = design
)

# -----------------------------
# 12. Fit contrasts + empirical Bayes
# -----------------------------

fit2 <- contrasts.fit(
  fit,
  contrasts
)

fit2 <- eBayes(fit2)

# -----------------------------
# 13. Save DE results
# -----------------------------

contrast_names <- colnames(contrasts)

for (contrast_name in contrast_names) {

  results <- topTable(
    fit2,
    coef = contrast_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )

  results$Probe_ID <- rownames(results)

  # Put probe ID first
  results <- results[, c(
    "Probe_ID",
    setdiff(colnames(results), "Probe_ID")
  )]

  output_file <- file.path(
    results_dir,
    paste0("DE_", contrast_name, ".csv")
  )

  write.csv(
    results,
    output_file,
    row.names = FALSE
  )

  cat(
    "\nSaved:",
    output_file,
    "\n"
  )
}

# -----------------------------
# 14. Summary
# -----------------------------

cat("\n========================================\n")
cat("Differential expression analysis complete\n")
cat("========================================\n")

cat("Samples:", nrow(expr_filtered), "\n")
cat("Probes tested:", ncol(expr_filtered), "\n")
cat("Contrasts:", length(contrast_names), "\n")

cat("\n===== SANITY CHECKS =====\n")

cat("\nSample counts:\n")
print(table(metadata$subtype))

cat("\nExpression dimensions:\n")
print(dim(expr))

cat("\nFiltered expression dimensions:\n")
print(dim(expr_filtered))

cat("\nMissing values:\n")
print(sum(is.na(expr_filtered)))

cat("\nDesign dimensions:\n")
print(dim(design))

cat("\nDesign matrix:\n")
print(design)

cat("\nContrast names:\n")
print(colnames(contrasts))