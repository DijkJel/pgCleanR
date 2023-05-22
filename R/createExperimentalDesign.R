#' Create experimental design data from a proteinGroups file
#'
#' Given a proteomics dataset \code{pg}, this function extracts information about the experimental design
#' and returns a data frame with columns for sample label, condition, and replicate. The sample labels are
#' taken from the LFQ.intensity columns of the dataset.
#'
#' @param pg A proteomics dataset.
#'
#' @return A data frame with columns for label, condition, and replicate.
#'
#' @examples
#' design <- createExperimentalDesign(pg)
#' @importFrom utils capture.output
#' @export
createExperimentalDesign = function(pg){

  lfq = grep('LFQ.intensity', names(pg))
  samples = gsub('LFQ.intensity.', '', names(pg)[lfq])
  condition = gsub('_.*$', '', samples)
  replicate = gsub('.*_', '', samples)
  
  df = data.frame(label = samples,
                  condition = condition,
                  replicate = replicate)
  
  df$replicate = as.vector(sapply(unique(df$condition), function(x){
    1:length(df$condition[df$condition == x])
  }))
  
  message(paste0(utils::capture.output(df), collapse = "\n"))
  return(df)
}
