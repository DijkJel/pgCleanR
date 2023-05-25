#' Get statistics for padj and ratio values after multiple rounds of imputation.
#'
#' Calculate the mean and standard deviation of p.adj values and ratios after multiple imputation.
#' In addition, it reports if proteins are significantly differentially expressed if both p.adj and ratio
#' meet the criteria set by padj_threshold and fc_threshold in more than fraction_cutoff * 100 % of the imputations
#'
#' @param matrix_list A list of matrices.
#' @param padj_threshold The p-value adjusted threshold (default: 0.05).
#' @param fc_threshold The fold change threshold (default: 2).
#' @param fraction_cutoff The fraction cutoff for significance (default: 0.8). In the default case, proteins need to be considered
#' significant in more than 80% of the rounds of imputations
#'
#' @return A list containing the following components:
#' \describe{
#' \item{padj_fractions}{A matrix of adjusted p-value fractions.}
#' \item{ratio_fractions}{A matrix of ratio fractions.}
#' \item{results}{A data frame with statistics for padj, ratio, and significance values}
#' }
#'
#'import matrixStats
#' @export
getImputationStats = function(matrix_list, padj_threshold = 0.05, fc_threshold = 2, fraction_cutoff = 0.8){
  
  m = matrix_list[[1]][[1]]
  
  padj_colnames = colnames(m)
  ratio_colnames = colnames(matrix_list[[2]][[1]])
  
  out_padj = vector(mode = 'list', length = nrow(m))
  names(out_padj) = rownames(m)
  out_ratio = vector(mode = 'list', length = nrow(m))
  names(out_ratio) = rownames(m)
  
  out_padj_fraction = matrix(, nrow = nrow(m), ncol = ncol(m))
  out_ratio_fraction = matrix(, nrow = nrow(m), ncol = ncol(m))
  
  rnames = rownames(m)
  
  for (p in 1:nrow(m)) {
    
    m = matrix_list[[1]][[1]]
    m_new_padj = matrix(, nrow = length(matrix_list[[1]]), ncol = ncol(m))
    m_new_ratio = matrix(, nrow = length(matrix_list[[1]]), ncol = ncol(m))
    for(m in 1:length(matrix_list[[1]])) {
      m_new_padj[m,] = matrix_list[[1]][[m]][p,]
      m_new_ratio[m,] = matrix_list[[2]][[m]][p,]
    }
    out_padj[[p]] = m_new_padj
    out_ratio[[p]] = m_new_ratio
    out_padj_fraction[p,] = apply(m_new_padj, 2, function(x) {sum(x >= padj_threshold) / nrow(m_new_padj)})
    out_ratio_fraction[p,] = apply(m_new_ratio, 2, function(x) {sum(abs(x) <= log2(fc_threshold)) / nrow(m_new_ratio)})
    
  }
  
  rownames(out_padj_fraction) = rnames
  rownames(out_ratio_fraction) = rnames

  if (ncol(out_padj[[1]]) == 1){
    mean_padj = as.matrix(sapply(out_padj, colMeans))
    sd_padj = as.matrix(sapply(out_padj, matrixStats::colSds))
    mean_ratio = as.matrix(sapply(out_ratio, colMeans))
    sd_ratio = as.matrix(sapply(out_ratio, matrixStats::colSds))
  }
  else{
    mean_padj = t(sapply(out_padj, colMeans))
    sd_padj = t(sapply(out_padj, matrixStats::colSds))
    mean_ratio = t(sapply(out_ratio, colMeans))
    sd_ratio = t(sapply(out_ratio, matrixStats::colSds))
  }
  

  
  colnames(mean_padj) = gsub('p.adj', 'mean_p.adj', padj_colnames)
  colnames(mean_ratio) = gsub('ratio', 'mean_ratio', ratio_colnames)
  colnames(sd_padj) = gsub('mean', 'sd', colnames(mean_padj))
  colnames(sd_ratio) = gsub('mean', 'sd', colnames(mean_ratio))
  
  zip_columns = function(a, b, c, d){
    
    l = list(a, b, c, d)
    df = as.data.frame(do.call(cbind, l))
    return(df)
  }
  
  results_stats = mapply(zip_columns, as.data.frame(mean_padj),
                         as.data.frame(sd_padj),
                         as.data.frame(mean_ratio),
                         as.data.frame(sd_ratio))
  
  results_stats = do.call(cbind, results_stats)
  
  
  colnames(results_stats) = mapply(function(a, b, c, d) c(a, b, c, d),
                                   colnames(mean_padj),
                                   colnames(sd_padj),
                                   colnames(mean_ratio),
                                   colnames(sd_ratio))
  
  
  #mean_results = cbind(mean_padj, mean_ratio)
  rownames(results_stats) = rownames(mean_padj)
  
  check_sig = function(x, y){ifelse(x < (1-fraction_cutoff) & y < (1-fraction_cutoff), T, F)}
  
  sig = mapply(check_sig, as.data.frame(out_padj_fraction), as.data.frame(out_ratio_fraction))
  colnames(sig) = gsub('p.adj', 'significant', padj_colnames)
  
  results = cbind(as.data.frame(results_stats), as.data.frame(sig))
  
  return(list(padj_fractions = out_padj_fraction,
              ratio_fractions = out_ratio_fraction,
              results = results))
  
  return(list(padj_fractions = out_padj_fraction,
              ratio_fractions = out_ratio_fraction,
              results_stats = results_stats,
              significant_hits = sig))
  
  #return(list(out_padj, out_ratio, out_padj_fraction, out_ratio_fraction, results_stats))
}
