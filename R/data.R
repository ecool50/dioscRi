#' Sample CyTOF Markers
#'
#' A character vector containing the names of 5 key CyTOF markers used
#' in the example datasets.
#'
#' @format A character vector of length 5 containing:
#' \describe{
#'   \item{CD3}{T cell marker}
#'   \item{CD4}{Helper T cell marker}
#'   \item{CD8}{Cytotoxic T cell marker}
#'   \item{CD19}{B cell marker}
#'   \item{CD14}{Monocyte marker}
#' }
#' @source Synthetic data
#' @examples
#' data(sample_markers)
#' print(sample_markers)
"sample_markers"

#' Sample CyTOF Single-Cell Data
#'
#' A dataset containing single-cell CyTOF measurements for 1,000 cells
#' from 10 samples with 5 markers.
#'
#' @format A data frame with 1,000 rows and 8 columns:
#' \describe{
#'   \item{cell_id}{Unique cell identifier (1-1000)}
#'   \item{sample_id}{Numeric sample identifier (1-10)}
#'   \item{batch}{Batch identifier (Batch1 or Batch2)}
#'   \item{cell_type}{Cell type assignment (CD4_T, CD8_T, B_cells, Monocytes, NK_cells)}
#'   \item{CD3}{CD3 marker expression (asinh transformed)}
#'   \item{CD4}{CD4 marker expression (asinh transformed)}
#'   \item{CD8}{CD8 marker expression (asinh transformed)}
#'   \item{CD19}{CD19 marker expression (asinh transformed)}
#'   \item{CD14}{CD14 marker expression (asinh transformed)}
#' }
#' @source Synthetic data generated with realistic cell type distributions
#' @examples
#' data(sample_cytof_data)
#' str(sample_cytof_data)
#' table(sample_cytof_data$cell_type)
#' table(sample_cytof_data$sample_id)
"sample_cytof_data"

#' Sample Clinical Data
#'
#' Clinical metadata for the 10 samples in the CyTOF dataset.
#'
#' @format A data frame with 10 rows and 4 columns:
#' \describe{
#'   \item{sample_id}{Numeric sample identifier (1-10)}
#'   \item{Outcome}{Binary outcome variable (factor: 0 or 1)}
#'   \item{Age}{Patient age in years}
#'   \item{Gender}{Patient gender (factor: M or F)}
#' }
#' @source Synthetic clinical data
#' @examples
#' data(sample_clinicaldata)
#' table(sample_clinicaldata$Outcome)
#' summary(sample_clinicaldata$Age)
"sample_clinicaldata"

#' Sample Logit-Transformed Cell Type Proportions
#'
#' A matrix of logit-transformed cell type proportions for each sample.
#' Used as features for group lasso modeling.
#'
#' @format A data frame with 10 rows (samples) and 5 columns (cell types):
#' \describe{
#'   \item{B_cells_logit}{Logit-transformed B cell proportion}
#'   \item{CD4_T_logit}{Logit-transformed CD4 T cell proportion}
#'   \item{CD8_T_logit}{Logit-transformed CD8 T cell proportion}
#'   \item{Monocytes_logit}{Logit-transformed monocyte proportion}
#'   \item{NK_cells_logit}{Logit-transformed NK cell proportion}
#' }
#' @details Proportions are logit-transformed using log(p/(1-p)) after
#' adjusting extreme values (0 → 0.001, 1 → 0.999).
#' @source Derived from sample_cytof_data
#' @examples
#' data(sample_data_logit)
#' dim(sample_data_logit)
#' head(sample_data_logit)
"sample_data_logit"

#' Sample Mean Marker Expressions
#'
#' A matrix of mean marker expressions per cell type for each sample.
#' Used as features for group lasso modeling.
#'
#' @format A data frame with 10 rows (samples) and 25 columns.
#' Column names follow the pattern: CellType_Marker
#' (e.g., "CD4_T_CD3", "B_cells_CD19").
#' @details Each value represents the mean expression of a specific marker
#' within a specific cell type for a given sample.
#' @source Derived from sample_cytof_data
#' @examples
#' data(sample_markerMeans)
#' dim(sample_markerMeans)
#' colnames(sample_markerMeans)[1:5]
"sample_markerMeans"

#' Sample Hierarchical Tree Structure
#'
#' A hierarchical tree structure generated from clustering cell type proportions.
#' Used for defining overlapping groups in group lasso.
#'
#' @format A list containing:
#' \describe{
#'   \item{tree}{A list with a data component containing tree structure}
#'   \item{tree$data}{Data frame with node information including clusters}
#'   \item{order}{Numeric vector with hierarchical clustering order}
#' }
#' @source Generated using hierarchical clustering on sample_data_logit
#' @examples
#' data(sample_tree)
#' names(sample_tree)
#' nrow(sample_tree$tree$data)
"sample_tree"

#' Sample Group Structure for Overlapping Group Lasso
#'
#' A list defining overlapping groups of variables for the group lasso penalty.
#' Each element contains indices of variables that belong to a group.
#'
#' @format A list of 10 numeric vectors, where each vector contains
#' variable indices (1-30) that form a group. Variables can appear in
#' multiple groups (overlapping).
#' @details Groups 1-3 contain proportion features (indices 1-5),
#' groups 4-8 contain cell-type-specific marker means,
#' and groups 9-10 contain marker-specific means across cell types.
#' @source Manually defined based on biological relationships
#' @examples
#' data(sample_groups)
#' length(sample_groups)
#' sapply(sample_groups, length)
"sample_groups"

#' Sample SingleCellExperiment Object
#'
#' A SingleCellExperiment object containing the CyTOF data with multiple
#' assays and metadata.
#'
#' @format A SingleCellExperiment object with:
#' \describe{
#'   \item{Dimensions}{5 markers × 1,000 cells}
#'   \item{Assays}{
#'     \itemize{
#'       \item norm: Asinh-transformed expressions
#'       \item raw: Original expressions (reverse-transformed)
#'       \item scaled: Z-score normalized expressions
#'     }
#'   }
#'   \item{colData}{Cell metadata including sample_id, batch, and cell_type}
#'   \item{rowData}{Marker metadata including marker names and classes}
#' }
#' @source Created from sample_cytof_data
#' @examples
#' data(sample_cytof_sce)
#' dim(sample_cytof_sce)
#' assayNames(sample_cytof_sce)
#' colData(sample_cytof_sce)[1:5,]
"sample_cytof_sce"