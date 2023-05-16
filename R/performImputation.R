#' Perform Imputation on Intensity Data
#'
#' This function takes log2-transformed and normalized LFQ intensity data as input and performs imputation to handle missing values.
#' Samples are split per conditions to determine whether MAR or MNAR is used. If all values for a condition are missing or in the lowest 10% of intensity values
#' MNAR is applied. Otherwise MAR is applied. MAR is performed by MLE, MNAR by minProb method.
#'
#' @param intensities numeric matrix of intensity data with samples in columns and features in rows.
#' @param expDesign A data frame specifying the experimental design with at least two columns: 'Sample.ID' and 'Condition'.
#' @return A numeric matrix with imputed intensity values for missing data.
#'
#' @import MsCoreUtils
#' @rawNamespace import(stats, except = c(start, end, smooth))
#' @export
performImputation = function(intensities, expDesign){

  conditions = unique(expDesign$condition)
  intensities = lapply(conditions, function(x){ints = intensities[,grep(x, colnames(intensities))]})

  imputed_intensities = lapply(intensities, function(x){
    
    
    q_10 = stats::quantile(as.vector(x), 0.1, na.rm = T)
    impute_type = apply(x, 1, function(y){
      if ((all(is.na(y)) | mean(y, na.rm = T) < q_10)){'MNAR'}
      else {'MAR'}
    })
    
    x = as.matrix(x)
    
    
    mnar = MsCoreUtils::impute_matrix(x, method = 'MinProb', q = 0.01)
    mar = MsCoreUtils::impute_mle(x)
    mnar = as.data.frame(t(mnar))
    mar = as.data.frame(t(mar))
    
    zip_impute = function(mnar, mar, mar_mnar){
      if (mar_mnar == 'MAR'){return(mar)}
      else {return(mnar)}
    }
    
    imputed_intensities = mapply(zip_impute, mnar, mar, impute_type)
  })
  
  intensities = do.call(rbind, imputed_intensities)
  intensities = t(as.data.frame(intensities))
  colnames(intensities) = expDesign$label
  
  
  return(intensities)
}