#' Perform Multiple Imputations
#'
#' Wrapper function that uses imputeMultiple and getImputationStats to
#' perform multiple imputations using the given parameters and returns the results of the imputations.
#'
#' @param pg A cleaned maxquant proteingroups data frame
#' @param expDesign An experimental design created with createExperimentalDesign, or with similar structure.
#' @param n The number of imputations to perform.
#' @param padj_threshold The p-value adjusted threshold for significance (default: 0.05).
#' @param fc_threshold The fold change threshold for significance (default: 2).
#' @param missing_threshold The missing data threshold for imputation (default: 1).
#' @param fraction_cutoff The fraction cutoff (default: 0.8).
#' 
#' @seealso \code{\link{performImputation}}, \code{\link{imputeMultiple}},
#'   \code{\link{getImputationStats}}
#'
#' @return The results of the imputation statistics.
#' @export
#'
#'

performMultipleImputations = function(pg, expDesign, n, padj_threshold = 0.05, fc_threshold = 2, missing_threshold = 1, fraction_cutoff = 0.8){
  
  multiple_imputations = imputeMultiple(pg, expDesign, n, padj_threshold, fc_threshold, missing_threshold)
  statistics = getImputationStats(multiple_imputations, padj_threshold, fc_threshold, fraction_cutoff)
  
  return(statistics$results)
}