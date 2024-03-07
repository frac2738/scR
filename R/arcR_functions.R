
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Subset an arcR projects based on sampleID value
#'
#' Given a sample name, the function returns an archR project containing only that 
#' sample.
#' @param arcR_project Directory in which the project will be located.  Default is "getwd()"
#' @param ident_value Name of the project
#' @export
#'
#' @return An archR project
#'
#' @concept ATAC-Seq
#'
#' @examples
#' \dontrun{
#' # If using built-in directory structure.
#' new_project <- archR_subset_by_sample(exp_Integrated,"C1")
#' }
#'
archR_subset_by_sample <- function(arcR_project, ident_value){
   
   idxSample <- BiocGenerics::which(arcR_project$Sample %in% ident_value)
   cellsSample <- arcR_project$cellNames[idxSample]
   tmp_dataset <- arcR_project[cellsSample, ]
   tmp_dataset <- addImputeWeights(tmp_dataset)
   
   return(tmp_dataset)
}
