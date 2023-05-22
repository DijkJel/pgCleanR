#' Impute Multiple
#'
#' Returns a list containing n matrices with p.adjust values and ratios for each pairwise comparison in experimentalDesign
#'
#' @param pg A cleaned proteinGroups file
#' @param expDesign The experimental design matrix as produced by createExperimentalDesign (or similar)
#' @param n The number of iterations for imputation.
#' @param padj_threshold The p-value adjustment threshold (default: 0.05).
#' @param fc_threshold The fold change threshold (default: 2).
#' @param missing_threshold The missing value threshold (default: 1).
#'
#' @return A list containing n matrices with p.adjust values and ratios for each comparison.
#'
#' @import DEP SummarizedExperiment
#' @export

imputeMultiple = function(pg, expDesign, n, padj_threshold = 0.05, fc_threshold = 2, missing_threshold = 1) {
  #Returns a list containing n matrices with p.adjust values  and ratios for each comparison made with add_rejections
  
  pg = DEP::make_unique(pg, "Gene.names", "Protein.IDs", delim = ";")
  data_se = DEP::make_se(pg, grep('LFQ', colnames(pg)), expDesign)
  data_filt = DEP::filter_missval(data_se, thr = missing_threshold)
  data_norm = DEP::normalize_vsn(data_filt)
  intensities = SummarizedExperiment::assay(data_norm)
  
  l_padj = vector(mode = 'list', length = n)
  l_ratio = vector(mode = 'list', length = n)
  
  for (i in 1:n) {
    
    #data_imputed = DEP::impute(dataset, fun = 'MinProb', q=0.01)
    data_imputed = performImputation(intensities, expDesign)
    SummarizedExperiment::assays(data_norm, withDimnames = F)[[1]] = data_imputed
    data_imputed = data_norm
    data_diff_all_contrasts <- DEP::test_diff(data_imputed, type = "all")
    dep <- DEP::add_rejections(data_diff_all_contrasts, alpha = padj_threshold, lfc = fc_threshold) 
    res = DEP::get_results(dep)
    r_names = res$name
    res_padj = as.matrix(res[,grep('p.adj', names(res)),drop = F])
    res_ratio = as.matrix(res[,grep('ratio', names(res)), drop = F])
    
    rownames(res_padj) = r_names
    rownames(res_ratio) = r_names
    res_padj = res_padj[order(rownames(res_padj)),,drop = F]
    res_ratio = res_ratio[order(rownames(res_ratio)),,drop = F]
    
    l_padj[[i]] = res_padj
    l_ratio[[i]] = res_ratio
  }
  
  return(list(l_padj, l_ratio))
  
}