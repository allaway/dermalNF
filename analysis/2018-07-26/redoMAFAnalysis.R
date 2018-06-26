### need to redo 

source("../../bin/WGSData.R")

allmafs <- getMAFs(mut.type='all')

#'getMutationSummary opens all maf files and loads into list
#'@param allmafs - list of MAFs to summarize
#'@return list of tables
#'
muts <- getMutationSummary(allmafs=allmafs)
