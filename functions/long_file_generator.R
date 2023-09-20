# Script to build custom long files using flat files
# Designed to use 'comorbidities.tsv' type input

# Depends
library(data.table)
library(stringr)

# Define long file generator
make_long_watchit <- function(def,cov_name,remove_missing_age=TRUE){
  out <- list(); n <- 1
  
  query <- def[cov==cov_name]$cov.code
  scope <- unique(def$cov.code.type)
  
  if ('MED' %in% scope){
    med_query <- def[c(cov==cov_name & cov.code.type=='MED')]$cov.code
    med_regex <- paste0(med_query,collapse='|')
  }
  
  for (i in 1:13){
    path <- list.files('/WATCH-IT_RPDR_files/RPDR_update2_2022')
    dir <- paste0('/WATCH-IT_RPDR_files/RPDR_update2_2022/',path[i],'/')
    dia_file <- list.files(dir)[str_detect(list.files(dir),'Dia')]
    dia <- fread(paste0(dir,dia_file),quote='',sep='|')
    output <- dia[Code %in% query,c('EMPI','Date','Code','Code_Type','Inpatient_Outpatient')]
    output[,Date := as.Date(Date,format='%m/%d/%Y')]
    if ('CPT' %in% scope){
      prc_file <- list.files(dir)[str_detect(list.files(dir),'Prc')]
      prc <- fread(paste0(dir,prc_file),quote='',sep='|')
      prc_output <- prc[Code %in% query,c('EMPI','Date','Code','Code_Type','Inpatient_Outpatient')]
      prc_output[,Date := as.Date(Date,format='%m/%d/%Y')]
      output <- rbind(output,prc_output)
    }
    if ('MED' %in% scope){
      med_file <- list.files(dir)[str_detect(list.files(dir),'Med')]
      med <- fread(paste0(dir,med_file),quote='',sep='|')
      med[,Medication := tolower(Medication)]
      setnames(med,c('Medication_Date','Medication'),c('Date','Code'))
      med_output <- med[str_detect(Code,regex(med_regex)),c('EMPI','Date','Code','Code_Type','Inpatient_Outpatient')]
      med_output[,Date := as.Date(Date,format='%m/%d/%Y')]
      output <- rbind(output,med_output)
    }
    out[[n]] <- output
    print(paste0("Just finished chunk: ",n))
    n <- n + 1
  }
  
  # Output processing
  ## Collapse list
  final <- do.call(rbind,out)
  ## Label the disease
  final[,cov := paste0(cov_name)]
  ## Set column order and row order (ascending per person)
  setorder(final,'cov')
  setkey(final,EMPI,Date)
  ## Remove rows with missing date (cannot use)
  if (remove_missing_age == TRUE){final <- final[!is.na(Date)]}
  ## Output
  return(final)
}