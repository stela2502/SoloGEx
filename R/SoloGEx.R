

#' @importFrom methods new
#' @importFrom stats sd
#' @importFrom utils write.table
#' @importFrom corrplot corrplot
#' @importFrom viridis viridis
#' @importFrom stats cor
#' @importFrom grDevices heat.colors

setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
#' SoloGEx S4 Class for Gene Expression Analysis
#'
#' The `SoloGEx` class is a minimal framework for storing gene expression data
#' along with sample metadata. It supports calculation of per-group statistics
#' (mean, SD, CV) and comparison of singleton (single-sample) datasets against
#' reference groups. Designed for exploratory analysis and hypothesis generation.
#'
#' @slot expr Numeric matrix of expression values. Rows correspond to genes,
#'   columns correspond to samples.
#' @slot colData Data frame of sample metadata. Each row corresponds to a column
#'   in \code{expr}. Must include columns identifying biological groups if
#'   computing group statistics.
#' @slot group_stats Optional matrix of per-gene, per-group statistics computed
#'   by \code{\link{compute_group_stats}}. Columns are of the form
#'   \code{<GroupName>_mean}, \code{<GroupName>_sd}, \code{<GroupName>_cv}.
#' @slot singleton_analysis Optional list containing results of
#'   \code{\link{analyze_singletons}}, including:
#'   \itemize{
#'     \item{\code{z_scores} – absolute Z-scores per gene vs closest group}
#'     \item{\code{closest_group} – reference group with minimum |Z| per gene}
#'     \item{\code{overall_dev} – mean minimum Z per sample}
#'   }
#'
#' @details
#' The `SoloGEx` class replaces an ExpressionSet-like object with a lightweight
#' list-based S4 implementation. Use the provided functions to:
#' \enumerate{
#'   \item{Create a new object: \code{SoloGEx(expr, colData)}}
#'   \item{Compute group statistics: \code{compute_group_stats()}}
#'   \item{Analyze singleton samples: \code{analyze_singletons()}}
#'   \item{Export combined results: \code{export_combined()} or \code{write_combined_file()}}
#' }
#'
#' @examples
#' # Create mock expression dataset
#' expr <- matrix(rnorm(12), nrow = 3, dimnames = list(paste0("Gene",1:3),
#'                                                     paste0("S",1:4)))
#' colData <- data.frame(SampleID = paste0("S",1:4),
#'                       BiologicalGroup = c("A","A","B","B"))
#' edat <- SoloGEx(expr, colData)
#'
#' # Compute group statistics
#' edat <- compute_group_stats(edat, group_col = "BiologicalGroup")
#'
#' @export
setClass(
  "SoloGEx",
  slots = list(
    expr = "matrix",                   # expression matrix genes x samples
    colData = "data.frame",            # sample metadata
    group_stats = "matrixOrNULL",      # group statistics (mean, SD, CV)
    singleton_analysis = "listOrNULL"  # singleton analysis results
    )
  )





#' SoloGEx Constructor
#'
#' Create a `SoloGEx` object to store gene expression data and sample metadata.
#'
#' @param expr A numeric matrix of expression values (genes × samples).
#' @param colData A data.frame of sample metadata; each row corresponds to a column in `expr`.
#'
#' @return A `SoloGEx` S4 object with slots:
#'   \itemize{
#'     \item{\code{expr} — expression matrix}
#'     \item{\code{colData} — sample metadata}
#'     \item{\code{group_stats} — initialized as NULL}
#'     \item{\code{singleton_analysis} — initialized as NULL}
#'   }
#'
#' @examples
#' # Example expression matrix
#' expr <- matrix(rnorm(12), nrow = 3, dimnames = list(paste0("Gene", 1:3),
#'                                                     paste0("S", 1:4)))
#' colData <- data.frame(SampleID = paste0("S", 1:4),
#'                       BiologicalGroup = c("A","A","B","B"))
#' edat <- SoloGEx(expr, colData)
#'
#' @export
SoloGEx <- function(expr, colData) {
  if (!is.matrix(expr)) expr <- as.matrix(expr)
  if (nrow(colData) != ncol(expr)) {
    stop(
      paste0(
        "ERROR in SoloGEx constructor: The number of rows in colData (",
        nrow(colData),
        ") must match the number of columns in expr (",
        ncol(expr),
        ").\n",
        "Each column of expr must correspond to a row in colData (sample metadata)."
        )
      )
  }  
  new("SoloGEx",
    expr = expr,
    colData = colData,
    group_stats = NULL,
    singleton_analysis = NULL)
}

setGeneric(
  "check_genes_match",
  function(obj1, obj2) standardGeneric("check_genes_match")
  )

#' Check gene consistency between two SoloGEx objects
#'
#' @param obj1 SoloGEx object (e.g., singletons)
#' @param obj2 SoloGEx object (e.g., larger reference set)
#' @return TRUE if genes match, stops with an informative error otherwise
#' @export
setGeneric(
  "check_genes_match",
  function(obj1, obj2) standardGeneric("check_genes_match")
  )

#' @describeIn check_genes_match Method for SoloGEx objects
setMethod(
  "check_genes_match",
  signature(obj1 = "SoloGEx", obj2 = "SoloGEx"),
  function( obj1, obj2){
    genes1 <- rownames(obj1@expr)
    genes2 <- rownames(obj2@expr)

    if (is.null(genes1) || is.null(genes2)) {
      stop("ERROR: One or both SoloGEx objects do not have rownames set for genes.")
    }

    if (!setequal(genes1, genes2)) {
      missing1 <- setdiff(genes1, genes2)
      missing2 <- setdiff(genes2, genes1)

      stop(
        "ERROR: Gene sets do not match between the two SoloGEx objects.\n",
        if(length(missing1) > 0) paste0("Genes in obj1 not in obj2: ", paste(missing1, collapse = ", "), "\n") else "",
        if(length(missing2) > 0) paste0("Genes in obj2 not in obj1: ", paste(missing2, collapse = ", "), "\n") else ""
        )
    }

    return(TRUE)
  }
)

setGeneric(
  "enforce_gene_order",
  function(target, reference) standardGeneric("enforce_gene_order")
  )


#' Reorder genes of target SoloGEx to match reference SoloGEx
#'
#' @param target SoloGEx object to reorder
#' @param reference SoloGEx object whose gene order will be used
#' @return target SoloGEx object with rows reordered to match reference
#' @export
setGeneric(
  "enforce_gene_order",
  function(target, reference) standardGeneric("enforce_gene_order")
  )

#' @describeIn enforce_gene_order Method for SoloGEx objects
setMethod(
  "enforce_gene_order",
  signature(target = "SoloGEx", reference = "SoloGEx"),
  function(target, reference ) {
  ref_genes <- rownames(reference@expr)
  target_genes <- rownames(target@expr)
  
  if (is.null(ref_genes) || is.null(target_genes)) {
    stop("ERROR: rownames must be set for both SoloGEx objects.")
  }
  
  if (!all(ref_genes %in% target_genes)) {
    stop("ERROR: Not all genes from reference exist in target.")
  }
  
  # Reorder rows
  target@expr <- target@expr[ref_genes, , drop = FALSE]
  
  return(target)
}
)


#' Compute per-gene, per-group statistics for a SoloGEx object
#'
#' This function calculates the mean, standard deviation (SD), and
#' coefficient of variation (CV = SD / mean) for each gene across
#' replicates within each biological group of a reference dataset.
#'
#' @param edat A \code{SoloGEx} S4 object containing the expression matrix
#'   in the \code{expr} slot and sample metadata in the \code{colData} slot.
#' @param group_col A character string specifying the column name in
#'   \code{colData} that identifies the biological groups.
#'
#' @return A \code{SoloGEx} object with the \code{group_stats} slot filled.
#'   The \code{group_stats} matrix has \code{genes x 3*groups} dimensions,
#'   with columns for each group: \code{<GroupName>_mean}, \code{<GroupName>_sd},
#'   \code{<GroupName>_cv}.
#'
#' @details
#' Each row corresponds to a gene.  
#' For each group:
#' \describe{
#'   \item{\code{<GroupName>_mean}}{ — mean expression across replicates}
#'   \item{\code{<GroupName>_sd}}{ — standard deviation across replicates}
#'   \item{\code{<GroupName>_cv}}{ — coefficient of variation (SD / mean)}
#' }
#' 
#' Interpretation:
#' \describe{
#'   \item{Low CV indicates a gene is stable across replicates in this group.}
#'   \item{High CV indicates a gene is noisy; such genes should be interpreted
#'         with caution in downstream analyses.}
#' }
#' 
#' @examples
#' # Create a mock SoloGEx object
#' expr <- matrix(rnorm(12), nrow = 3, dimnames = list(paste0("Gene",1:3),
#'                                                     paste0("S",1:4)))
#' colData <- data.frame(SampleID = paste0("S",1:4),
#'                       BiologicalGroup = c("A","A","B","B"))
#' edat <- SoloGEx(expr, colData)
#' group_stats <- compute_group_stats(edat, group_col = "BiologicalGroup")
#' head(group_stats)
#'
#' @export
setGeneric("compute_group_stats", function(edat, group_col) {
  standardGeneric("compute_group_stats")
  })

#' @describeIn compute_group_stats Method for SoloGEx objects
setMethod(
  "compute_group_stats",
  signature(edat = "SoloGEx"),
  function(edat, group_col) {

    # Check expr is non-empty
    if (nrow(edat@expr) == 0 || ncol(edat@expr) == 0) {
      stop("ERROR: Expression matrix is empty.")
    }
    
    # Check colData exists and has the specified column
    if (is.null(edat@colData) || !(group_col %in% colnames(edat@colData))) {
      stop(
        "ERROR: group_col '", group_col, 
        "' not found in colData. Available columns: ", 
        paste(colnames(edat@colData), collapse = ", ")
        )
    }
    
    # Check group column is non-empty
    if (any(is.na(edat@colData[[group_col]]))) {
      warning(
        "WARNING: NA values found in group_col '", group_col, 
        "'. These samples will be ignored."
        )
    }
    
    # Check that number of columns in expr matches rows in colData
    if (ncol(edat@expr) != nrow(edat@colData)) {
      stop(
        "ERROR: Number of columns in expr (", ncol(edat@expr), 
        ") does not match number of rows in colData (", nrow(edat@colData), ")."
        )
    }

    expr <- edat@expr
    groups <- edat@colData[[group_col]]
    ugs <- unique(groups)

    n_genes <- nrow(expr)
    n_groups <- length(ugs)

  # Initialize result matrix
  stats <- matrix(NA, n_genes, n_groups * 3)
  colnames(stats) <- as.vector(rbind(
    paste0(ugs, "_mean"),
    paste0(ugs, "_sd"),
    paste0(ugs, "_cv")
    ))
  rownames(stats) <- rownames(expr)
  
  for (i in seq_along(ugs)) {
    g <- ugs[i]
    idx <- which(groups == g)
    sub <- expr[, idx, drop = FALSE]
    
    mean_vals <- rowMeans(sub, na.rm = TRUE)
    sd_vals <- apply(sub, 1, sd, na.rm = TRUE)
    cv_vals <- sd_vals / mean_vals
    
    stats[, (i-1)*3 + 1] <- mean_vals
    stats[, (i-1)*3 + 2] <- sd_vals
    stats[, (i-1)*3 + 3] <- cv_vals
  }
  return(stats)
})



#' Analyze singleton samples against reference group statistics
#'
#' This function compares singleton (1-sample) datasets to a reference
#' dataset with group statistics. For each gene in each singleton sample:
#' it computes absolute Z-scores relative to each reference group,
#' identifies the "closest" group (smallest |Z|), and calculates an
#' overall deviation per sample.
#'
#' @param single_edat A \code{SoloGEx} object containing singleton samples.
#'   Expression data should be in \code{@expr}, sample metadata in \code{@colData}.
#' @param group_edat A \code{SoloGEx} object containing the reference dataset.
#'   If group statistics are not yet computed, they will be calculated using
#'   \code{compute_group_stats()} internally.
#' @param group_col Character string specifying the column in \code{colData} 
#'   that defines biological groups.
#'
#' @return A \code{SoloGEx} object with the \code{singleton_analysis} slot filled:
#'   \describe{
#'     \item{\code{z_scores}}{Absolute Z-score of each gene vs closest group.}
#'     \item{\code{closest_group}}{Name of the reference group with minimum |Z| per gene.}
#'     \item{\code{overall_dev}}{Mean minimum Z-score per sample, indicating global deviation.}
#'   }
#'
#' @details
#' Steps performed internally:
#' \enumerate{
#'   \item{Check that genes in \code{single_edat} and \code{group_edat} match.}
#'   \item{Compute group statistics for \code{group_edat} if missing.}
#'   \item{For each gene in each singleton sample:}
#'     \itemize{
#'       \item{Compute Z-score vs mean of each reference group.}
#'       \item{Identify group with minimum absolute Z-score.}
#'     }
#'   \item{Compute overall deviation per sample as the mean of min-Z across genes.}
#' }
#'
#' Interpretation for lab scientists:
#' \itemize{
#'   \item{High |Z| values indicate genes unusually high or low relative to reference groups.}
#'   \item{A high \code{overall_dev} for a sample suggests it is globally different from all reference groups.}
#'   \item{Closest group per gene can guide hypotheses about which baseline the singleton most resembles.}
#' }
#'
#' @examples
#' # Create mock reference dataset
#' ref_expr <- matrix(rnorm(12), nrow=3, dimnames=list(paste0("Gene",1:3), paste0("R",1:4)))
#' ref_colData <- data.frame(SampleID = paste0("R",1:4),
#'                           BiologicalGroup = c("A","A","B","B"))
#' ref_edat <- SoloGEx(ref_expr, ref_colData)
#'
#' # Create singleton dataset
#' singleton_expr <- matrix(rnorm(6), nrow=3, dimnames=list(paste0("Gene",1:3), paste0("S",1:2)))
#' singleton_colData <- data.frame(SampleID = paste0("S",1:2))
#' single_edat <- SoloGEx(singleton_expr, singleton_colData)
#'
#' # Analyze
#' single_edat <- analyze_singletons(single_edat, ref_edat, group_col="BiologicalGroup")
#' single_edat@singleton_analysis
#'
#' @export
setGeneric("analyze_singletons", function(single_edat, group_edat, group_col) {
  standardGeneric("analyze_singletons")
  })


#' @describeIn analyze_singletons Method for SoloGEx objects
setMethod(
  "analyze_singletons",
  signature(single_edat = "SoloGEx", group_edat= "SoloGEx"), 
  function(single_edat, group_edat, group_col = NULL) {
    single_expr <- single_edat@expr
    n_genes <- nrow(single_expr)
    n_samples <- ncol(single_expr)

    check_genes_match(single_edat, group_edat )


    group_stats = compute_group_stats (group_edat, group_col )

    single_edat@group_stats = group_stats

    group_names <- unique(sub("_(mean|sd|cv)$", "", colnames(group_stats)[seq(1, ncol(group_stats), 3)]))


    z_scores <- matrix(NA, n_genes, n_samples, dimnames = list(rownames(single_expr), colnames(single_expr)))
    closest_group <- matrix("", n_genes, n_samples, dimnames = list(rownames(single_expr), colnames(single_expr)))

  # Loop through each singleton sample
  for (s in 1:n_samples) {
    min_z <- rep(Inf, n_genes)
    min_group <- rep(NA_character_, n_genes)
    
    # Compare to each reference group
    for (i in seq_along(group_names)) {
      g <- group_names[i]
      
      # Extract mean and SD for this group
      mean_vals <- group_stats[, (i-1)*3 + 1]
      sd_vals   <- group_stats[, (i-1)*3 + 2]
      
      # Replace zero SDs with NA to avoid division by zero
      sd_vals[sd_vals == 0] <- NA_real_
      
      # Compute z-scores
      z_tmp <- abs((single_expr[, s] - mean_vals) / sd_vals)
      
      # Update minimum Z and closest group safely
      update <- z_tmp < min_z
      update[is.na(update)] <- FALSE  # skip positions where z_tmp is NA
      
      min_z[update] <- z_tmp[update]
      min_group[update] <- g
    }
    
    z_scores[, s] <- min_z
    closest_group[, s] <- min_group
  }
  
  overall_dev <- colMeans(z_scores, na.rm = TRUE)
  
  single_edat@singleton_analysis <- list(
    z_scores = z_scores,
    closest_group = closest_group,
    overall_dev = overall_dev
    )
  
  return(single_edat)
})



#' Export combined table from a SoloGEx object
#'
#' Combines the expression matrix, group statistics, and singleton analysis
#' from a \code{SoloGEx} object into a single \code{data.frame} suitable for
#' downstream export (e.g., CSV) and interpretation by lab scientists.
#'
#' @param edat A \code{SoloGEx} S4 object. The object can have:
#'   \describe{
#'     \item{\code{@expr}}{ — expression matrix (genes x samples)}
#'     \item{\code{@group_stats}}{ — optional, per-group statistics (mean, SD, CV)}
#'     \item{\code{@singleton_analysis}}{ — optional, singleton analysis
#'           (Z-scores, closest group, overall deviation)}
#'   }
#'
#' @return A \code{data.frame} combining:
#'   \itemize{
#'     \item{Expression values (genes x samples)}
#'     \item{Group statistics (mean, SD, CV) if available}
#'     \item{Singleton analysis results (Z-scores, closest group, overall deviation) if available}
#'   }
#'
#' @details
#' The resulting table contains all genes (rows) and the following columns:
#' \enumerate{
#'   \item{Expression values for each sample}
#'   \item{Group statistics for each reference group (columns named <Group>_mean, <Group>_sd, <Group>_cv)}
#'   \item{Singleton analysis results:
#'     \describe{
#'       \item{Z-scores}](columns starting with "Z_")}
#'       \item{\code{closest_group}}{indicating the most similar reference group per gene}
#'       \item{\code{overall_dev}}{summarizing mean deviation per singleton sample}
#'     } }
#' }
#'
#' Interpretation for lab scientists:
#' \itemize{
#'   \item{High Z-scores indicate genes unusually different from reference groups.}
#'   \item{\code{closest_group} helps to hypothesize which baseline a singleton resembles.}
#'   \item{\code{overall_dev} indicates global deviation of a singleton sample.}
#' }
#'
#' @examples
#' # Assume 'single_edat' has singleton analysis and group stats computed
#' \dontrun{
#' combined_table <- export_combined(single_edat)
#' head(combined_table)
#' }
#'
#' @export
setGeneric("export_combined", function( edat ) {
  standardGeneric("export_combined")
  })

#' @describeIn export_combined Method for SoloGEx objects
setMethod(
  "export_combined",
  signature(edat = "SoloGEx"),
  function(edat) {
    # Start with expression matrix
    df <- as.data.frame(edat@expr)

    fc = compute_fc_pairwise( edat )
    df = cbind ( df, fc )

    # Append group statistics if available
    if (!is.null(edat@group_stats)) {
      df <- cbind(df, as.data.frame(edat@group_stats))
    }

    
    
    # Append singleton analysis if available
    if (!is.null(edat@singleton_analysis)) {
      # Z-scores
      z_df <- as.data.frame(edat@singleton_analysis$z_scores)
      colnames(z_df) <- paste0("Zscore_", colnames(edat@singleton_analysis$z_scores))

      # Closest group
      closest_group <- as.data.frame(edat@singleton_analysis$closest_group)
      colnames(closest_group) <- paste0("closest_", colnames(edat@singleton_analysis$closest_group))

      # Overall deviation
      #overall_dev <- data.frame(overall_dev = edat@singleton_analysis$overall_dev)
      
      df <- cbind(df, z_df, closest_group)
    }
    
    return(df)
  }
  )


#' Write combined SoloGEx results to a file
#'
#' Exports the combined results from a \code{SoloGEx} object to a file.
#' This includes expression values, reference group statistics, and singleton
#' analysis results. The function is a wrapper around \code{export_combined()}
#' and standard \code{write.table()}.
#'
#' @param singleton A \code{SoloGEx} object with \code{singleton_analysis} slot filled.
#'   Expression data should be in \code{@expr}, sample metadata in \code{@colData}.
#' @param file Character string specifying the file path to write the table.
#' @param ... Additional arguments passed to \code{write.table()}, e.g.,
#'   \code{sep = ","}, \code{quote = FALSE}, \code{col.names = TRUE}.
#'
#' @return Invisibly returns \code{NULL}. The results are written to the specified file.
#'
#' @details
#' Before writing, the function checks that the \code{singleton_analysis} slot exists.
#' If missing, it stops with an informative error. The resulting file contains:
#' \itemize{
#'   \item{Expression values for each gene and sample}
#'   \item{Group statistics (mean, SD, CV) if available}
#'   \item{Singleton analysis results (Z-scores, closest group, overall deviation)}
#' }
#'
#' @examples
#' # Assume 'single_edat' has singleton analysis computed
#' \dontrun{
#' write_combined_file(single_edat, "singleton_results.csv", sep = ",", quote = FALSE)
#' }
#'
#' @export
setGeneric(
  "write_combined_file",
  function(singleton, file, ...) standardGeneric("write_combined_file")
  )

#' @describeIn write_combined_file Method for SoloGEx objects
setMethod(
  "write_combined_file",
  signature(singleton = "SoloGEx"),
  function(singleton, file, ...) {
    if (is.null(singleton@singleton_analysis) ) {
      stop("No singleton analysis found. Run analyze_singletons() first.")
    }
    
    combined <- export_combined(singleton)
    combined = cbind(rownames( combined), combined )
    colnames(combined)[1] = "GeneNames"
    write.table(combined, file = file, ..., row.names = FALSE)
  }
  )


# Internal helper: compute all pairwise log2 fold changes for singleton samples
compute_fc_pairwise <- function(single_edat) {

  single_expr <- single_edat@expr
  n_genes <- nrow(single_expr)
  n_samples <- ncol(single_expr)
  sample_names <- colnames(single_expr)
  
  fc_list <- list()
  
  for (i in seq_len(n_samples - 1)) {
    for (j in seq((i + 1), n_samples)) {
      fc_name <- paste0("logFC_", sample_names[i], "_vs_", sample_names[j])
      
      # Extract the two columns
      x <- single_expr[, i]
      y <- single_expr[, j]
      
      # Initialize fold change vector with NAs
      fc_vec <- rep(NA_real_, n_genes)
      
      # Compute log2 fold change only where valid (non-NA and denominator != 0)
      valid <- !is.na(x) & !is.na(y) & y != 0
      fc_vec[valid] <- log2(x[valid] / y[valid])
      
      # Assign to list
      fc_list[[fc_name]] <- fc_vec
    }
  }
  
  # Combine into a data.frame with genes as rows
  fc_df <- as.data.frame(fc_list, row.names = rownames(single_expr))
  return(fc_df)
}


#' Compute correlation of singleton samples to reference samples (S4 method)
#'
#' Calculates correlation between each singleton sample and all reference samples.
#'
#' @param singleton SoloGEx object containing singleton expression data
#' @param reference SoloGEx object containing reference expression data
#' @param method Correlation method: "pearson" or "spearman" (default: "pearson")
#' 
#' @return SoloGEx object with correlations stored in @singleton_analysis$correlation_to_reference
#' @export
setGeneric(
  "compute_corr_to_reference",
  function(singleton, reference, method = "pearson") standardGeneric("compute_corr_to_reference")
)

#' @describeIn compute_corr_to_reference Method for SoloGEx objects
setMethod(
  "compute_corr_to_reference",
  signature(singleton = "SoloGEx", reference = "SoloGEx"),
  function(singleton, reference, method = "pearson") {

    # Check that singleton and reference have expression matrices
    if (is.null(singleton@expr) || ncol(singleton@expr) == 0 || nrow(singleton@expr) == 0) {
      stop("ERROR: Singleton SoloGEx object has no expression data.")
    }
    if (is.null(reference@expr) || ncol(reference@expr) == 0 || nrow(reference@expr) == 0) {
      stop("ERROR: Reference SoloGEx object has no expression data.")
    }

    # Ensure singleton_analysis exists
    if (is.null(singleton@singleton_analysis)) {
      stop("ERROR: singleton_analysis slot is empty. Run analyze_singletons() first.")
    }

    # Ensure gene sets match
    check_genes_match(singleton, reference)
    singleton <- enforce_gene_order(singleton, reference)
    
    # Compute correlation matrix
    single_expr <- singleton@expr
    ref_expr <- reference@expr
    cor_mat <- matrix(NA, nrow = ncol(single_expr), ncol = ncol(ref_expr),
                      dimnames = list(colnames(single_expr), colnames(ref_expr)))
    
    for (s in seq_len(ncol(single_expr))) {
      for (r in seq_len(ncol(ref_expr))) {
        cor_mat[s, r] <- cor(single_expr[, s], ref_expr[, r], method = method, use = "pairwise.complete.obs")
      }
    }
    
    # Store correlations in singleton_analysis
    singleton@singleton_analysis$correlation_to_reference <- cor_mat
    return(singleton)
  }
)

#' Plot correlation of singleton samples to reference samples using corrplot
#'
#' @param singleton SoloGEx object with @singleton_analysis$correlation_to_reference
#' @param method Plot method for corrplot: "circle", "color", or "number" (default: "circle")
#' @param title Plot title (default: "Singleton vs Reference Correlation")
#' @param ... Additional arguments passed to corrplot
#' @export
setGeneric(
  "plot_corr_heatmap",
  function(singleton, method = "circle", title = "Singleton vs Reference Correlation", ...) 
    standardGeneric("plot_corr_heatmap")
)

#' @describeIn plot_corr_heatmap Method for SoloGEx objects
setMethod(
  "plot_corr_heatmap",
  signature(singleton = "SoloGEx"),
  function(singleton, method = "circle", title = "Singleton vs Reference Correlation", ...) {
    # Check that correlation matrix exists
    if (is.null(singleton@singleton_analysis$correlation_to_reference)) {
      stop("ERROR: correlation_to_reference not found. Run compute_corr_to_reference() first.")
    }
    
    cor_mat <- singleton@singleton_analysis$correlation_to_reference
    
    # Ensure corrplot package is available
    if (!requireNamespace("corrplot", quietly = TRUE)) {
      stop("The 'corrplot' package is required. Install with install.packages('corrplot').")
    }
    
    # Set color palette: high correlation dark, low correlation light
    if (!requireNamespace("viridis", quietly = TRUE)) {
      colors <- heat.colors(100)
    } else {
      colors <- rev(viridis::viridis(100))
    }
    
    # Create the plot
    corrplot::corrplot(cor_mat,
                       method = method,
                       col = colors,
                       tl.col = "black",
                       addCoef.col = if (method == "number") "white" else NULL,
                       number.cex = 0.8,
                       mar = c(0,0,1,0),
                       title = title,
                       ...)
  }
)

################################################################################
# USAGE EXAMPLE
# (Replace 'ref_expr', 'ref_info', 'single_expr', 'single_info' with your actual data)
################################################################################

# # Reference dataset
# ref_edat <- SoloGEx(expr = ref_expr, colData = ref_info)
# ref_edat <- compute_group_stats(ref_edat, "BiologicalGroup")
# 
# # Singleton dataset
# single_edat <- SoloGEx(expr = single_expr, colData = single_info)
# single_edat <- analyze_singletons(single_edat, ref_edat@group_stats)
# 
# # Export full combined table for lab scientists
# \dontrun{
# combined_table <- export_combined(single_edat)
# head(combined_table)
# }
################################################################################

################################################################################
# EXCESSIVE EXPLANATION FOR LAB SCIENTISTS
################################################################################
# 1) SoloGEx:
#    - Stores your expression data and metadata in one object.
#    - Ensures that expression columns match metadata rows.
#    - Slots for computed group stats and singleton analysis.
#
# 2) compute_group_stats:
#    - Calculates baseline metrics from your reference dataset.
#    - Mean, SD, and CV per gene per group.
#    - Allows us to understand which genes are stable vs noisy.
#    - This is critical because single-sample comparisons rely on these variances.
#
# 3) analyze_singletons:
#    - Compares singletons to all reference groups.
#    - Z-scores tell you how extreme a gene’s expression is relative to the closest reference group.
#    - Closest group per gene informs about similarity to known biological groups.
#    - overall_dev tells you how globally different a sample is.
#    - Interpretation:
#        * High Z → gene potentially interesting, may need validation.
#        * Low Z → gene behaves as expected.
#        * High overall_dev → sample deviates strongly from all reference groups.
#
# 4) export_combined:
#    - Combines everything in a single, easy-to-read table.
#    - Useful for downstream reporting, plotting, or sharing with wet-lab colleagues.
#    - Can be imported into Excel, R, Python, or other analysis tools.
#
# IMPORTANT:
#    - This is exploratory/pilot-level analysis. No inferential statistics (p-values) are computed.
#    - Single-sample results are hypothesis-generating only.
#    - Always validate high-Z genes or unusual samples with replicates if possible.
################################################################################
