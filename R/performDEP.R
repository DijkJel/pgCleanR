#' Perform Differential Expression Analysis using DEP 
#'
#' This function performs differential expression analysis on protein expression data and prioritizes significant genes/proteins based on adjusted p-values and log fold change (lfc).
#' @param pg A cleaned maxquant proteinGroups data frame containing LFQ protein expression data with gene names and protein IDs.
#' @param expDesign A data frame specifying the experimental design with at least two columns: 'Sample.ID' and 'Condition'.
#' @param thr A numeric value specifying the threshold for filtering out missing values (default is 1).
#' @return A data frame containing gene/protein IDs and other identifiers with significant genes/proteins prioritized based on adjusted p-values and lfc.
#' @examples
#' \dontrun{performDEP(pg, expDesign)}
#' @import DEP SummarizedExperiment
#' @importFrom assertthat assert_that
#' @export

performDEP = function(pg, expDesign, thr = 1, padj_cutoff = 0.05, fc_cutoff = 2){
  
  
  stopifnot(checkExperimentalDesign(expDesign, pg))
  
  pg = DEP::make_unique(pg, "Gene.names", "Protein.IDs", delim = ";")
  data_se = DEP::make_se(pg, grep('LFQ', colnames(pg)), expDesign)
  data_filt = DEP::filter_missval(data_se, thr = thr)
  data_norm = DEP::normalize_vsn(data_filt)
  intensities = SummarizedExperiment::assay(data_norm)
  
  intensities_imputed = performImputation(intensities, expDesign = expDesign)
  SummarizedExperiment::assays(data_norm, withDimnames = F)[[1]] = intensities_imputed
  data_imp = data_norm
  
  data_diff_all_contrasts = DEP::test_diff(data_imp, type = "all")
  dep = DEP::add_rejections(data_diff_all_contrasts, alpha = padj_cutoff, lfc = log2(fc_cutoff)) 
  results = DEP::get_results(dep)
  
  ids = pg[,grep('Protein.ID|Ensembl|Entrez|RefSeq', colnames(pg))]
  ids = merge(ids, results, all.x = F, all.y = T, by.x = 'Protein.IDs', by.y = 'ID')
  ids = ids[order(ids$significant, decreasing = T),]
  
  out = list(results = ids, intensities = intensities_imputed)
  return(out)
}
