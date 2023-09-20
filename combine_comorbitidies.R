# Script to combine/widen disease long files for watch-it

# Depends
library(data.table)
library(plyr)
library(stringr)

# Load master
full <- fread(file='watchit_full_101922.csv')

###### AF
af_long <- fread(file='af_long.csv')
af_prc <- af_long[c(Code_Type=='CPT' | Code=='99.61'),.SD[1],by=EMPI]
af_dia_in <- af_long[c((Code_Type == 'ICD9' | Code_Type == "ICD10") & Inpatient_Outpatient=="Inpatient"),
                     .SD[1],by=EMPI]
af_any_code <- af_long[,.SD[2],by=EMPI][!is.na(Code)]
af_any_unique <- unique(af_long[EMPI %in% af_any_code$EMPI],by='EMPI')$EMPI
af_any_code <- af_long[EMPI %in% af_any_unique,.SD[1],by=EMPI]

## Add definitions in sequence
setkey(full,empi); setkey(af_prc,EMPI); setkey(af_dia_in,EMPI); setkey(af_any_code,EMPI)

full <- full[af_any_code,af_dia_date := i.Date]
full <- full[af_prc,af_prc_date := i.Date]
full <- full[af_dia_in,af_dia_in_date := i.Date]

## Min is AF date
full[,af_date := pmin(af_dia_date,af_prc_date,af_dia_in_date,na.rm=TRUE)]

## Other arrhythmias
oa_long <- fread('oa_long.csv')
oa_long <- oa_long[,.SD[1],by=EMPI][!is.na(Code)]
setkey(oa_long,EMPI)
full <- full[oa_long,oa_date := i.Date]

## CAD
cad_long <- fread('cad_long.csv')
cad_long <- cad_long[,.SD[1],by=EMPI][!is.na(Code)]
setkey(cad_long,EMPI)
full <- full[cad_long,cad_date := i.Date]

## PAD
pad_long <- fread('pad_long.csv')
pad_long <- pad_long[,.SD[1],by=EMPI][!is.na(Code)]
setkey(pad_long,EMPI)
full <- full[pad_long,pad_date := i.Date]

## HF
hf_long <- fread('hf_long.csv')
hf_long <- hf_long[Inpatient_Outpatient=='Inpatient',.SD[1],by=EMPI][!is.na(Code)]
setkey(hf_long,EMPI)
full <- full[hf_long,hf_date := i.Date]

## CVA
cva_long <- fread('cva_long.csv')
cva_long <- cva_long[,.SD[1],by=EMPI][!is.na(Code)]
setkey(cva_long,EMPI)
full <- full[cva_long,cva_date := i.Date]

## Dementia
dementia_long <- fread('dementia_long.csv')
dementia_long <- dementia_long[,.SD[1],by=EMPI][!is.na(Code)]
setkey(dementia_long,EMPI)
full <- full[dementia_long,dementia_date := i.Date]

## DM
dm_long <- fread('dm_long.csv')
dm_long <- dm_long[,.SD[1],by=EMPI][!is.na(Code)]
setkey(dm_long,EMPI)
full <- full[dm_long,dm_date := i.Date]

## HTN
htn_long <- fread('htn_long.csv')
htn_long <- htn_long[,.SD[1],by=EMPI][!is.na(Code)]
setkey(htn_long,EMPI)
full <- full[htn_long,htn_date := i.Date]

## HLD
hld_long <- fread('hld_long.csv')
hld_long <- hld_long[,.SD[1],by=EMPI][!is.na(Code)]
setkey(hld_long,EMPI)
full <- full[hld_long,hld_date := i.Date]

# Obesity
obesity_long <- fread('obesity_long.csv')
obesity_long <- obesity_long[,.SD[1],by=EMPI][!is.na(Code)]
setkey(obesity_long,EMPI)
full <- full[obesity_long,obesity_date := i.Date]

## Pulmonary
pulm_long <- fread('pulm_long.csv')
pulm_long <- pulm_long[,.SD[1],by=EMPI][!is.na(Code)]
setkey(pulm_long,EMPI)
full <- full[pulm_long,pulm_date := i.Date]

## Malignancy
malignancy_long <- fread('malignancy_long.csv')
malignancy_long <- malignancy_long[,.SD[1],by=EMPI][!is.na(Code)]
setkey(malignancy_long,EMPI)
full <- full[malignancy_long,malignancy_date := i.Date]

############ VARIABLE CLEANUP
# Dates to ages
full[,':='(af_age = (af_date-birth_date)/365.25,
           oa_age = (oa_date-birth_date)/365.25,
           cad_age = (cad_date-birth_date)/365.25,
           pad_age = (pad_date-birth_date)/365.25,
           hf_age = (hf_date-birth_date)/365.25,
           cva_age = (cva_date-birth_date)/365.25,
           dementia_age = (dementia_date-birth_date)/365.25,
           dm_age = (dm_date-birth_date)/365.25,
           htn_age = (htn_date-birth_date)/365.25,
           hld_age = (hld_date-birth_date)/365.25,
           obesity_age = (obesity_date-birth_date)/365.25,
           pulm_age = (pulm_date-birth_date)/365.25,
           malignancy_age = (malignancy_date-birth_date)/365.25)]

# Study start date
full[,study_start := as.Date('2021-11-16',format='%Y-%m-%d')]
full[,':='(study_start_age = (as.numeric(study_start)-as.numeric(birth_date))/365.25)]

# Prevalent disease markers
full[,':='(prev_af = ifelse(c(!is.na(af_age) & (af_age <= study_start_age)),1,0),
           prev_oa = ifelse(c(!is.na(oa_age) & (oa_age <= study_start_age)),1,0),
           prev_cad = ifelse(c(!is.na(cad_age) & (cad_age <= study_start_age)),1,0),
           prev_pad = ifelse(c(!is.na(pad_age) & (pad_age <= study_start_age)),1,0),
           prev_hf = ifelse(c(!is.na(hf_age) & (hf_age <= study_start_age)),1,0),
           prev_cva = ifelse(c(!is.na(cva_age) & (cva_age <= study_start_age)),1,0),
           prev_dementia = ifelse(c(!is.na(dementia_age) & (dementia_age <= study_start_age)),1,0),
           prev_dm = ifelse(c(!is.na(dm_age) & (dm_age <= study_start_age)),1,0),
           prev_htn = ifelse(c(!is.na(htn_age) & (htn_age <= study_start_age)),1,0),
           prev_hld = ifelse(c(!is.na(hld_age) & (hld_age <= study_start_age)),1,0),
           prev_obesity = ifelse(c(!is.na(obesity_age) & (obesity_age <= study_start_age)),1,0),
           prev_pulm = ifelse(c(!is.na(pulm_age) & (pulm_age <= study_start_age)),1,0),
           prev_malignancy = ifelse(c(!is.na(malignancy_age) & (malignancy_age <= study_start_age)),1,0))]

# Write out full
write.csv(full,file='full_phenos.csv',row.names=F)
