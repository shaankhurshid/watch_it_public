# Script to generate comorbidities/baseline chars for WATCH-IT

# Depends
library(data.table)
library(stringr)
library(plyr)

#### STAGE 1: CLEANUP

# Load existing datasets
full <- fread(file='fullcohort.csv')
ecare <- fread(file='ecare.csv')
redcap <- fread(file='redcap6_29.csv')

# Formatting fixes for ecare dataset
names(ecare) <- c('V1','patient','empi','dob','portal_status','ok_to_contact','date_contacted')
ecare[,':='(last_name = str_remove_all(tolower(str_extract(patient,regex('^\\w+',ignore_case = T))),'\\s+'),
            first_name = str_remove_all(tolower(str_extract(str_extract(patient,regex('\\, \\w+',ignore_case = T)),'\\w+')),'\\s+'))]

# Formatting fixes for redcap dataset
redcap[,dob_formatted := as.Date(paste0(str_extract(dob,'\\d+\\/\\d+\\/'),'19',str_extract(dob,'\\d+$')),format='%m/%d/%Y')]
redcap[,dob_formatted := ifelse(as.numeric(dob_formatted,origin='1970-01-01') <= -21915,
                                as.numeric(dob_formatted,origin='1970-01-01') + 36525,dob_formatted)] #Born before 1910
redcap[,dob_formatted := as.Date(dob_formatted,origin='1970-01-01')]

redcap[,':='(last_name = str_remove_all(tolower(last_name),'\\s+'),
             first_name = str_remove_all(tolower(first_name),'\\s+'),
             last_partial = str_remove_all(tolower(substr(last_name,1,3)),'\\s+'))]
redcap <- redcap[!c(is.na(dob_formatted) | is.na(last_name) | last_name=='' | last_name=='NA' | last_name=='na' |
                      is.na(first_name) | first_name=='' | first_name=='NA' | first_name=='na')] # 14986 - 1498 = 13488

# Select complete entries only
redcap_complete <- redcap[survey_complete=='Complete']
redcap_complete_unique <- unique(redcap_complete,by=c('first_name','last_name','dob_formatted'))

redcap_incomplete <- redcap[survey_complete=='Incomplete']
redcap_incomplete_unique <- unique(redcap_incomplete,by=c('first_name','last_name','dob_formatted'))

redcap <- rbind(redcap_complete,redcap_incomplete)
redcap <- unique(redcap,by=c('first_name','last_name','dob_formatted'))

# Merge full and ecare
setkey(full,empi); setkey(ecare,empi)
full[ecare,':='(in_ecare = 1, portal_status = i.portal_status,
                ok_to_contact = i.ok_to_contact, date_contacted = i.date_contacted,
                last_name = i.last_name, first_name = i.first_name,
                last_partial = str_remove_all(tolower(substr(last_name,1,3)),'\\s+'))]

# Merge with redcap
setkey(full,birth_date,last_partial); setkey(redcap,dob_formatted,last_partial)
full[redcap,':='(in_redcap = 1, redcap_id = i.id, survey_complete = i.survey_complete,
                 redcap_last = i.last_name, redcap_first = i.first_name)]

# In redcap failed to merge
not_found = redcap[!(id %in% full$redcap_id)]
write.csv(not_found,file='redcap_not_merged.csv',row.names=F) #1728
write.csv(full,file='merged_10142022.csv',row.names=F)

#### STAGE 2: PHENOTYPES
## Use linkers to get study-specific linker_ids
c3po_linker <- fread(file='patient_linker_file_c3po062120.csv')
ewoc_linker <- fread(file='c3po_cardiology_id_linker.txt')

## Merge
setkey(c3po_linker,EMPI); setkey(ewoc_linker,EMPI); setkey(full,empi)
full[c3po_linker,c3po_linker_id := i.linker_id]
full[ewoc_linker,ewoc_linker_id := i.linker_id]

## Create long files for Dems
#### C3PO
out <- list()
for (i in 1:13){
  path <- list.files('WATCH-IT_RPDR_files/RPDR_update2_2022')
  dir <- paste0('WATCH-IT_RPDR_files/RPDR_update2_2022/',path[i],'/')
  dem_file <- list.files(dir)[str_detect(list.files(dir),'Dem')]
  out[[i]] <- fread(paste0(dir,dem_file),quote='',sep='|')
  watch_dem <- do.call(rbind,out)
  print(paste0("Just finished chunk: ",i))
}
write.csv(watch_dem,file='watchit_dem_long.csv',row.names=F)
setkey(watch_dem,EMPI); setkey(full,empi)
full[watch_dem, ':=' (race = i.Race_Group, ethnicity = i.Ethnic_Group)]
full[,race_ethnicity := ifelse(ethnicity=='HISPANIC','hispanic',
                               ifelse(race=='White','white',
                                      ifelse(race=='Black','black',
                                             ifelse(c(race=='Asian' | race=='Native Hawaiian or Other Pacific Islander'),'asian',
                                                    ifelse(race=='Two or More','mixed',
                                                           ifelse(c(race=='Unknown/Missing' | is.na(race) | race=='Declined'),NA,'other'))))))]

# Write out full
write.csv(full,file='watchit_full_101922.csv',row.names=F)

#### STAGE 3A: CREATE LONG FILES FOR DISEASES
# Load comorbs
comorbs <- fread(file='comorbidities.tsv')
af_procs <- fread(file='AF_Procedures_mod.txt')

## AF
# Isolate to AF and get it into the right format
af_procs[,':='(cov = 'af_aflt')]
setnames(af_procs,c('codeType','code'),c('cov.code.type','cov.code'))
af_procs <- af_procs[,c('cov','cov.code','cov.code.type')]

af <- rbind(af_procs,comorbs[cov=='af_aflt',c('cov','cov.code','cov.code.type')])
af <- unique(af,by='cov.code')

# Longs
af_long <- make_long_watchit(def=af,cov_name='af_aflt')
write.csv(af_long,file='af_long.csv',row.names=F)

## Other arrhythmias
other_arrhythmias <- comorbs[cov=='othArrhyth',c('cov','cov.code','cov.code.type')]

# Longs
oa_long <- make_long_watchit(def=other_arrhythmias,cov_name='othArrhyth')
write.csv(oa_long,file='oa_long.csv',row.names=F)

## CAD
cad <- comorbs[cov=='cad',c('cov','cov.code','cov.code.type')]

# Longs
cad_long <- make_long_watchit(def=cad,cov_name='cad')
write.csv(cad_long,file='cad_long.csv',row.names=F)

## PAD
pad <- comorbs[cov=='pad',c('cov','cov.code','cov.code.type')]

# Longs
pad_long<- make_long_watchit(def=pad,cov_name='pad')
write.csv(pad_long,file='pad_long.csv',row.names=F)

## HF
hf <- comorbs[cov=='heartFailure',c('cov','cov.code','cov.code.type')]

# Longs
hf_long<- make_long_watchit(def=hf,cov_name='heartFailure')
write.csv(hf_long,file='hf_long.csv',row.names=F)

## CVA
cva <- comorbs[cov=='cvaTia',c('cov','cov.code','cov.code.type')]

# Longs
cva_long <- make_long_watchit(def=cva,cov_name='cvaTia')
write.csv(cva_long,file='cva_long.csv',row.names=F)

## Dementia
dementia <- comorbs[cov=='dementia',c('cov','cov.code','cov.code.type')]

# Longs
dementia_long <- make_long_watchit(def=dementia,cov_name='dementia')
write.csv(dementia_long,file='dementia_long.csv',row.names=F)

## DM
dm <- comorbs[cov=='dm',c('cov','cov.code','cov.code.type')]

# Longs
dm_long <- make_long_watchit(def=dm,cov_name='dm')
write.csv(dm_long,file='dm_long.csv',row.names=F)

## HTN
htn <- comorbs[cov=='htn',c('cov','cov.code','cov.code.type')]

# Longs
htn_long <- make_long_watchit(def=htn,cov_name='htn')
write.csv(htn_long,file='htn_long.csv',row.names=F)

## HLD
hld <- comorbs[cov=='lipid',c('cov','cov.code','cov.code.type')]

# Longs
hld_long <- make_long_watchit(def=hld,cov_name='lipid')
write.csv(hld_long,file='hld_long.csv',row.names=F)

## Obesity
obesity <- comorbs[cov=='obesity',c('cov','cov.code','cov.code.type')]

# Longs
obesity_long <- make_long_watchit(def=obesity,cov_name='obesity')
write.csv(obesity_long,file='obesity_long.csv',row.names=F)

## Pulmonary
pulm <- comorbs[cov=='pulmonary',c('cov','cov.code','cov.code.type')]

# Longs
pulm_long <- make_long_watchit(def=pulm,cov_name='pulmonary')
write.csv(pulm_long,file='pulm_long.csv',row.names=F)

## Malignancy
malignancy <- comorbs[cov=='malignancy',c('cov','cov.code','cov.code.type')]

# Longs
malignancy_long <- make_long_watchit(def=malignancy,cov_name='malignancy')
write.csv(malignancy_long,file='malignancy_long.csv',row.names=F)




