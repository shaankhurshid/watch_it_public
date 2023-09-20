# Script to analyze WATCH-IT data

# Depends
library(data.table)
library(plyr)
library(stringr)
library(sktools)

# Load full
full <- fread(file='full_phenos.csv')
redcap <- fread(file='redcap6_29.csv')
contacted <- full[!is.na(date_contacted) | in_redcap==1]
non_respond <- full[!is.na(date_contacted)][is.na(in_redcap)]

# Remove columns that will result in duplication
full <- full[,.SD,.SDcols=names(full)[!(names(full) %in% c('survey_complete','in_redcap','redcap_last','redcap_first'))]]

# Merge
setkey(full,redcap_id); setkey(redcap,id)
redcap_merged <- full[redcap,nomatch=0]

# Add sex and zip code
watchit_dem_long <- fread('watchit_dem_long.csv')
setkey(watchit_dem_long,EMPI); setkey(redcap_merged,empi); setkey(non_respond,empi)
redcap_merged[watchit_dem_long,sex := i.Gender_Legal_Sex]
non_respond[watchit_dem_long,sex := i.Gender_Legal_Sex]
redcap_merged[watchit_dem_long,zip := i.Zip_code]
non_respond[watchit_dem_long,zip := i.Zip_code]

# Add income
income <- fread(file='2020_5y_zip_income.csv')
names(income) <- paste0(income[1])
income <- income[2:nrow(income)]
names(income)[names(income)=='Estimate!!Households!!Median income (dollars)'] <- 'median_income'
names(income)[names(income)=='Geographic Area Name'] <- 'zip_code'
income$zip_code <- as.numeric(str_remove(income$zip_code,'ZCTA5 '))
income[,median_income := as.numeric(median_income)]
income <- income[!(zip_code==0)]

redcap_merged[,zip_numeric := as.numeric(zip)]
non_respond[,zip_numeric := as.numeric(zip)]
setkey(redcap_merged,zip_numeric); setkey(non_respond,zip_numeric); setkey(income,zip_code)
redcap_merged[income,med_income := i.median_income]
non_respond[income,med_income := i.median_income]

# Analyze
setnames(redcap_merged,'Do you currently use a wearable device (e.g., Smartwatch or fitness tracker)?','user')

# User and non-user
user <- redcap_merged[!is.na(user) & user=='Yes']
non_user <- redcap_merged[!is.na(user) & user=='No']

# Write outs
write.csv(redcap_merged,file='redcap_merged_110222.csv',row.names=F)
write.csv(user,file='user_110222.csv',row.names=F)
write.csv(non_user,file='non_user_110222.csv',row.names=F)
write.csv(non_respond,file='non_respond_110222.csv',row.names=F)

# Model
#### A. User vs. Non user
# Shallow copy
redcap_merged_analysis <- redcap_merged

# Imputations
redcap_merged_analysis[is.na(med_income)]$med_income <- median(redcap_merged_analysis$med_income,na.rm=T)

# Variable formatting
redcap_merged_analysis[,':='(cards_cat = factor(ifelse(EWOC==1,'cardiology','non_cardiology'),levels=c('non_cardiology','cardiology')),
                             race_binary = ifelse(!is.na(race) & race=="White",'white','non_white'),
                             ethnicity_binary = ifelse(!is.na(ethnicity) & ethnicity=='HISPANIC','hispanic','non-hispanic'),
                             age10 = study_start_age/10,
                             income_cat = quantilize(med_income,4),
                             user_binary = ifelse(user=='Yes',1,0))]

# Models
use <- glm(user_binary ~ age10 + sex + race_binary + income_cat + cards_cat + prev_af + 
             prev_cad + prev_dm + prev_hld + prev_htn + prev_hf + prev_oa + prev_obesity + prev_pad + prev_pulm + prev_cva,
           data=redcap_merged_analysis,family='binomial')
or <- data.table(names=c('intercept','age10','sex','race_binary','income_cat','cards_cat','prev_af','prev_cad',
                         'prev_dm','prev_hld','prev_htn','prev_hf','prev_oa','prev_obesity','prev_pad','prev_pulm','prev_cva'),
                 or=exp(use$coefficients),lower=exp(confint(use)[,1]),upper=exp(confint(use)[,2]))
or_plot <- or[!(names=='intercept')]

# Plots
# Plot
color <- ifelse(c(or_plot$lower > 1 | or_plot$upper < 1),'#cb181d','darkgray')
pdf('user.pdf',height=4,width=6,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,15,1,1))
plot(x=or_plot$or,y=seq(15.5,0.5,-1),xlim=c(0.5,1.5),ylim=c(0,16),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=19,col=color,cex=2,bty='n')
axis(1,cex.axis=1.6,at=seq(0.5,1.5,0.1),pos=0)
axis(2,at=seq(15.5,0.5,-1),labels=FALSE,cex=1.8)
mtext('Factor',side=2,line=14,cex=2)
mtext('Adjusted odds ratio',side=1,line=3,cex=2)
segments(or_plot$lower,seq(15.5,0.5,-1),or_plot$upper,seq(15.5,0.5,-1),col=color,lwd=2.2)
segments(1,0,1,15,col='black',lty=5)

plot_names <- c('Age (per 10y increase)','Male sex','White race','Income (per quartile)','Cardiology clinic',
                'Atrial fibrillation','Coronary disease','Diabetes','Hyperlipidemia',
                'Hypertension','Heart failure','Other arrhythmias','Obesity',
                'Peripheral artery disease','Pulmonary disease','Stroke')
for (i in seq(15.5,0.5,-1)){
  n <- 16-i
  mtext(paste0(plot_names[i+0.5]),
        side=2,line=1,las=2,cex=1.5,at=n,
        font=ifelse(color=='darkgray',NA,2)[i+0.5])
}

dev.off()

#### B. Potential vs no potential
# Shallow copy
potential_analysis <- redcap_merged_analysis[c(is.na(user_binary) | user_binary==0)]

# Variable formatting
setnames(potential_analysis,'Would you consider using a wearable device if one were provided to you at no cost?','potential_user')
potential_analysis[,':='(potential_user_binary = ifelse(is.na(potential_user),NA,
                                                        ifelse(potential_user=='No',0,1)))]

# Models
pot_use <- glm(potential_user_binary ~ age10 + sex + race_binary + income_cat + cards_cat + prev_af + 
                 prev_cad + prev_dm + prev_hld + prev_htn + prev_hf + prev_oa + prev_obesity + prev_pad + prev_pulm + prev_cva,
               data=potential_analysis,family='binomial')
or <- data.table(names=c('intercept','age10','sex','race_binary','income_cat','cards_cat','prev_af','prev_cad',
                         'prev_dm','prev_hld','prev_htn','prev_hf','prev_oa','prev_obesity','prev_pad','prev_pulm','prev_cva'),
                 or=exp(pot_use$coefficients),lower=exp(confint(pot_use)[,1]),upper=exp(confint(pot_use)[,2]))
or_plot <- or[!(names=='intercept')]

# Plots
# Plot
color <- ifelse(c(or_plot$lower > 1 | or_plot$upper < 1),'#cb181d','darkgray')
pdf('pot_user.pdf',height=4,width=6,
    pointsize=3)
par(oma=c(1,1,1,1))
par(mar=c(4,15,1,1))
plot(x=or_plot$or,y=seq(15.5,0.5,-1),xlim=c(0,3),ylim=c(0,16),
     xaxt='n',yaxt='n',xlab='',ylab='',pch=19,col=color,cex=2,bty='n')
axis(1,cex.axis=1.6,at=seq(0,3,0.5),pos=0)
axis(2,at=seq(15.5,0.5,-1),labels=FALSE,cex=1.8)
mtext('Factor',side=2,line=14,cex=2)
mtext('Adjusted odds ratio',side=1,line=3,cex=2)
segments(or_plot$lower,seq(15.5,0.5,-1),or_plot$upper,seq(15.5,0.5,-1),col=color,lwd=2.2)
segments(1,0,1,16,col='black',lty=5)

plot_names <- c('Age (per 10y increase)','Male sex','White race','Income (per quartile)','Cardiology clinic',
                'Atrial fibrillation','Coronary disease','Diabetes','Hyperlipidemia',
                'Hypertension','Heart failure','Other arrhythmias','Obesity',
                'Peripheral artery disease','Pulmonary disease','Stroke')
for (i in seq(15.5,0.5,-1)){
  n <- 16-i
  mtext(paste0(plot_names[i+0.5]),
        side=2,line=1,las=2,cex=1.5,at=n,
        font=ifelse(color=='darkgray',NA,2)[i+0.5])
}

dev.off()

# Bar graph 1:Use
#count(redcap_merged_analysis$user_binary)
pdf(file='barplot_use.pdf',height=6,width=5,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(5,7,2,1))

coords <- barplot(c(6206,4915,183),
                  col=c('#4292c6','#ef6548','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,6500))

axis(2,at=seq(0,6500,500),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=3250,cex=2)
mtext('User',1,line=1.5,at=coords[1],cex=2)
mtext('Non User',1,line=1.5,at=coords[2],cex=2)
mtext('No Response',1,line=1.5,at=coords[3],cex=2)
dev.off()

# Bar graph 2: MD recommend
pdf(file='md_recommend.pdf',height=6,width=5,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(5,7,2,1))

coords <- barplot(c(867,532,9147),
                  col=c('#4292c6','darkgray','#ef6548'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,10000))

axis(2,at=seq(0,10000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=5000,cex=2)
mtext('Yes',1,line=1.5,at=coords[1],cex=2)
mtext('Unsure',1,line=1.5,at=coords[2],cex=2)
mtext('No',1,line=1.5,at=coords[3],cex=2)
dev.off()

# Bar graph 3: Telehealth
pdf(file='tele.pdf',height=6,width=5,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(9245,79,1258,722),
                  col=c('#4292c6','#a6bddb','#ef6548','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,10000))

axis(2,at=seq(0,10000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=5000,cex=2)
mtext('Yes',1,line=1.5,at=coords[1],cex=2)
mtext('Unsure',1,line=1.5,at=coords[2],cex=2)
mtext('No',1,line=1.5,at=coords[3],cex=2)
mtext('No Response',1,line=1.5,at=coords[4],cex=2)
dev.off()

# Bar graph 4: Brand
pdf(file='rand.pdf',height=6,width=5,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(5,7,2,1))

coords <- barplot(c(3485,1753,419,241,813),
                  col=c('#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,4000))

axis(2,at=seq(0,4000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=2000,cex=2)
mtext('Apple',1,line=1.5,at=coords[1],cex=2)
mtext('Fitbit',1,line=1.5,at=coords[2],cex=2)
mtext('Garmin',1,line=1.5,at=coords[3],cex=2)
mtext('Samsung',1,line=1.5,at=coords[4],cex=2)
mtext('Other',1,line=1.5,at=coords[5],cex=2)
dev.off()

# Bar graph 5: How many days (N=14 non-response)
#count(redcap_merged$`In the last 7 days, how many days did you wear your device for at least 12 hours?`)
pdf(file='days.pdf',height=6,width=5,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(5,7,2,1))

coords <- barplot(c(268,374,1084,4466),
                  col=c('#cb181d','#fc9272','#a6bddb','#3690c0'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,5000))

axis(2,at=seq(0,5000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=2500,cex=2)
mtext('0 days',1,line=1.5,at=coords[1],cex=2)
mtext('1-3 days',1,line=1.5,at=coords[2],cex=2)
mtext('4-6 days',1,line=1.5,at=coords[3],cex=2)
mtext('7 days',1,line=1.5,at=coords[4],cex=2)
dev.off()

# Bar graph 5A: Important (actual, N=15 non-response)
pdf(file='important_actual.pdf',height=6,width=5,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(5,7,2,1))

coords <- barplot(c(966,2075,1229,1921),
                  col=c('#cb181d','#a6bddb','#3690c0','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,3000))

axis(2,at=seq(0,3000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=1500,cex=2)
mtext('Not Important',1,line=1.5,at=coords[1],cex=2)
mtext('Moderate',1,line=1.5,at=coords[2],cex=2)
mtext('Very',1,line=1.5,at=coords[3],cex=2)
mtext('Unsure',1,line=1.5,at=coords[4],cex=2)
dev.off()

# Bar graph 5B: Important (potential, N=0 non-response)
pdf(file='important_potential.pdf',height=6,width=5,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(5,7,2,1))

coords <- barplot(c(228,1533,2328,597),
                  col=c('#cb181d','#a6bddb','#3690c0','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,3000))

axis(2,at=seq(0,3000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=1500,cex=2)
mtext('Not Important',1,line=1.5,at=coords[1],cex=2)
mtext('Moderate',1,line=1.5,at=coords[2],cex=2)
mtext('Very',1,line=1.5,at=coords[3],cex=2)
mtext('Unsure',1,line=1.5,at=coords[4],cex=2)
dev.off()

# Bar graph 6A: Next step (actual)
pdf(file='important_share_actual.pdf',height=6,width=9,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(5,7,2,1))

coords <- barplot(c(3995,1601,1096,5566,3329,3103,3058,1720,110,289),
                  col=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#e31a1c',
                        '#fb9a99','#fdbf6f','#ff7f00','#cab2d6','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,6000))

axis(2,at=seq(0,6000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=3000,cex=2)
mtext('Activity',1,line=1.5,at=coords[1],cex=2)
mtext('Anxiety',1,line=1.5,at=coords[2],cex=2)
mtext('Diet',1,line=1.5,at=coords[3],cex=2)
mtext('Heart rate',1,line=1.5,at=coords[4],cex=2)
mtext('Health',1,line=1.5,at=coords[5],cex=2)
mtext('Conditions',1,line=3,at=coords[5],cex=2)
mtext('Pulse',1,line=1.5,at=coords[6],cex=2)
mtext('Oximetry',1,line=3,at=coords[6],cex=2)
mtext('Sleep',1,line=1.5,at=coords[7],cex=2)
mtext('Weight',1,line=1.5,at=coords[8],cex=2)
mtext('Other',1,line=1.5,at=coords[9],cex=2)
mtext('None',1,line=1.5,at=coords[10],cex=2)
dev.off()

# Bar graph 6B: Next step (potential)
pdf(file='important_share_potential.pdf',height=6,width=9,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(5,7,2,1))

coords <- barplot(c(3434,2270,1915,4321,3626,3303,3385,2496,243,27),
                  col=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#e31a1c',
                        '#fb9a99','#fdbf6f','#ff7f00','#cab2d6','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,6000))

axis(2,at=seq(0,6000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=3000,cex=2)
mtext('Activity',1,line=1.5,at=coords[1],cex=2)
mtext('Anxiety',1,line=1.5,at=coords[2],cex=2)
mtext('Diet',1,line=1.5,at=coords[3],cex=2)
mtext('Heart rate',1,line=1.5,at=coords[4],cex=2)
mtext('Health',1,line=1.5,at=coords[5],cex=2)
mtext('Conditions',1,line=3,at=coords[5],cex=2)
mtext('Pulse',1,line=1.5,at=coords[6],cex=2)
mtext('Oximetry',1,line=3,at=coords[6],cex=2)
mtext('Sleep',1,line=1.5,at=coords[7],cex=2)
mtext('Weight',1,line=1.5,at=coords[8],cex=2)
mtext('Other',1,line=1.5,at=coords[9],cex=2)
mtext('None',1,line=1.5,at=coords[10],cex=2)
dev.off()

# Bar graph 6: No cost
#count(non_user$`Would you consider using a wearable device if one were provided to you at no cost?`)
pdf(file='free.pdf',height=6,width=7,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(5,7,2,1))

coords <- barplot(c(4336,350,146,266),
                  col=c('#3690c0','#a6bddb','#cb181d','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,5000))

axis(2,at=seq(0,4500,500),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=2250,cex=2)
mtext('Yes',1,line=1.5,at=coords[1],cex=2)
mtext('Yes if recommended',1,line=1.5,at=coords[2],cex=2)
mtext('by provider',1,line=3,at=coords[2],cex=2)
mtext('No',1,line=1.5,at=coords[3],cex=2)
mtext('No Response',1,line=1.5,at=coords[4],cex=2)
dev.off()

# Bar graph 7: Would share with doctor
pdf(file='would_share.pdf',height=6,width=5,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(4567,116,4,228),
                  col=c('#4292c6','#a6bddb','#ef6548','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,5000))

axis(2,at=seq(0,5000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=2500,cex=2)
mtext('Yes',1,line=1.5,at=coords[1],cex=2)
mtext('Unsure',1,line=1.5,at=coords[2],cex=2)
mtext('No',1,line=1.5,at=coords[3],cex=2)
mtext('No Response',1,line=1.5,at=coords[4],cex=2)
dev.off()

# Bar graph 8: Reassured
pdf(file='would_share.pdf',height=6,width=5,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(4567,116,4,228),
                  col=c('#4292c6','#a6bddb','#ef6548','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,5000))

axis(2,at=seq(0,5000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=2500,cex=2)
mtext('Yes',1,line=1.5,at=coords[1],cex=2)
mtext('Unsure',1,line=1.5,at=coords[2],cex=2)
mtext('No',1,line=1.5,at=coords[3],cex=2)
mtext('No Response',1,line=1.5,at=coords[4],cex=2)
dev.off()

# Bar graph 10: Doctor recommended wearable
pdf(file='md_recommend_wearable.pdf',height=6,width=5,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(868,536,9178,722),
                  col=c('#4292c6','#a6bddb','#ef6548','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,10000))

axis(2,at=seq(0,10000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=5000,cex=2)
mtext('Yes',1,line=1.5,at=coords[1],cex=2)
mtext('Unsure',1,line=1.5,at=coords[2],cex=2)
mtext('No',1,line=1.5,at=coords[3],cex=2)
mtext('No Response',1,line=1.5,at=coords[4],cex=2)
dev.off()

# Bar graph 11: Would you share with researchers?
pdf(file='would_share_research.pdf',height=6,width=6,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(3824,1835,67,214,266),
                  col=c('#4292c6','#a6bddb','#ef6548','darkgray','lightgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,4000))

axis(2,at=seq(0,4000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=2000,cex=2)
mtext('Extremely',1,line=1.5,at=coords[1],cex=2)
mtext('Somewhat',1,line=1.5,at=coords[2],cex=2)
mtext('Not at all',1,line=1.5,at=coords[3],cex=2)
mtext('Unsure',1,line=1.5,at=coords[4],cex=2)
mtext('No Response',1,line=1.5,at=coords[5],cex=2)
dev.off()

# Bar graph 12: Why wouldnt you use?
pdf(file='why_not.pdf',height=5,width=8,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(6,5,5,1))

coords <- barplot(c(41,43,15,51,57),
                  col=c('#1f78b4','#33a02c','#e31a1c',
                        '#ff7f00','#cab2d6'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,60))

axis(2,at=seq(0,60,10),las=2,cex.axis=1.8)
mtext('Count',2,line=4,at=30,cex=2)
mtext('Privacy',1,line=1.5,at=coords[1],cex=2)
mtext('Uncomfortable',1,line=1.5,at=coords[2],cex=2)
mtext('Would not know',1,line=1.5,at=coords[3],cex=2)
mtext('how to use it',1,line=3,at=coords[3],cex=2)
mtext('Not interested',1,line=1.5,at=coords[4],cex=2)
mtext('Other',1,line=1.5,at=coords[5],cex=2)
dev.off()

# Bar graph 13: Peace of mind (actual)
pdf(file='piece_of_mind_actual.pdf',height=5,width=8,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(3828,1984,104,70,220),
                  col=c('#4292c6','#a6bddb','#fb9a99','#ef6548','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,4000))

axis(2,at=seq(0,4000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=2000,cex=2)
mtext('Strongly Agree',1,line=1.5,at=coords[1],cex=2)
mtext('Agree',1,line=1.5,at=coords[2],cex=2)
mtext('Disagree',1,line=1.5,at=coords[3],cex=2)
mtext('Strongly Disagree',1,line=1.5,at=coords[4],cex=2)
mtext('No Response',1,line=1.5,at=coords[5],cex=2)
dev.off()

# Bar graph 13: Peace of mind (potential)
pdf(file='piece_of_mind_potential.pdf',height=5,width=8,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(2846,1594,116,359),
                  col=c('#4292c6','#a6bddb','#ef6548','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,4000))

axis(2,at=seq(0,4000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=2000,cex=2)
mtext('Strongly Agree',1,line=1.5,at=coords[1],cex=2)
mtext('Agree',1,line=1.5,at=coords[2],cex=2)
mtext('Disagree',1,line=1.5,at=coords[3],cex=2)
mtext('No Response',1,line=1.5,at=coords[4],cex=2)
dev.off()

# Bar graph 14: Peace of mind (potential)
pdf(file='piece_of_mind_potential.pdf',height=5,width=8,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(2846,1594,116,359),
                  col=c('#4292c6','#a6bddb','#ef6548','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,4000))

axis(2,at=seq(0,4000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=2000,cex=2)
mtext('Strongly Agree',1,line=1.5,at=coords[1],cex=2)
mtext('Agree',1,line=1.5,at=coords[2],cex=2)
mtext('Disagree',1,line=1.5,at=coords[3],cex=2)
mtext('No Response',1,line=1.5,at=coords[4],cex=2)
dev.off()

# Bar graph 15A: Anxiety (actual)
pdf(file='anxiety_actual.pdf',height=5,width=8,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(200,347,2052,3386,221),
                  col=c('#4292c6','#a6bddb','#fb9a99','#ef6548','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,4000))

axis(2,at=seq(0,4000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=2000,cex=2)
mtext('Strongly Agree',1,line=1.5,at=coords[1],cex=2)
mtext('Agree',1,line=1.5,at=coords[2],cex=2)
mtext('Disagree',1,line=1.5,at=coords[3],cex=2)
mtext('Strongly Disagree',1,line=1.5,at=coords[4],cex=2)
mtext('No Response',1,line=1.5,at=coords[5],cex=2)
dev.off()


# Bar graph 15B: Anxiety (potential)
pdf(file='anxiety_potential.pdf',height=5,width=8,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(219,419,3917,360),
                  col=c('#4292c6','#a6bddb','#ef6548','darkgray'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,4000))

axis(2,at=seq(0,4000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=2000,cex=2)
mtext('Strongly Agree',1,line=1.5,at=coords[1],cex=2)
mtext('Agree',1,line=1.5,at=coords[2],cex=2)
mtext('Disagree',1,line=1.5,at=coords[3],cex=2)
mtext('No Response',1,line=1.5,at=coords[4],cex=2)
dev.off()

# Bar graph 16A: Next step (actual)
pdf(file='next_step_actual.pdf',height=5,width=8,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(5444,1009,989,2865,2695,3903,170),
                  col=c('#1f78b4','#33a02c','#e31a1c',
                        '#fb9a99','#fdbf6f','#ff7f00','#cab2d6'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,6000))

axis(2,at=seq(0,6000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=3000,cex=2)
mtext('Contact',1,line=1.5,at=coords[1],cex=2)
mtext('provider',1,line=3,at=coords[1],cex=2)
mtext('Urgent care',1,line=1.5,at=coords[2],cex=2)
mtext('Telehealth',1,line=1.5,at=coords[3],cex=2)
mtext('Search',1,line=1.5,at=coords[4],cex=2)
mtext('internet',1,line=3,at=coords[4],cex=2)
mtext('Share results',1,line=1.5,at=coords[5],cex=2)
mtext('with family',1,line=3,at=coords[5],cex=2)
mtext('Track',1,line=1.5,at=coords[6],cex=2)
mtext('Symptoms',1,line=3,at=coords[6],cex=2)
mtext('Other',1,line=1.5,at=coords[7],cex=2)
dev.off()

# Bar graph 16B: Next step (potential)
pdf(file='next_step_potential.pdf',height=5,width=8,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(4315,791,704,2085,2051,2748,162),
                  col=c('#1f78b4','#33a02c','#e31a1c',
                        '#fb9a99','#fdbf6f','#ff7f00','#cab2d6'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,5000))

axis(2,at=seq(0,5000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=2500,cex=2)
mtext('Contact',1,line=1.5,at=coords[1],cex=2)
mtext('provider',1,line=3,at=coords[1],cex=2)
mtext('Urgent care',1,line=1.5,at=coords[2],cex=2)
mtext('Telehealth',1,line=1.5,at=coords[3],cex=2)
mtext('Search',1,line=1.5,at=coords[4],cex=2)
mtext('internet',1,line=3,at=coords[4],cex=2)
mtext('Share results',1,line=1.5,at=coords[5],cex=2)
mtext('with family',1,line=3,at=coords[5],cex=2)
mtext('Track',1,line=1.5,at=coords[6],cex=2)
mtext('Symptoms',1,line=3,at=coords[6],cex=2)
mtext('Other',1,line=1.5,at=coords[7],cex=2)
dev.off()

# Bar graph 17: MD communication
pdf(file='md_communication.pdf',height=5,width=8,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(2735,9398,5040,584,150),
                  col=c('#1f78b4','#33a02c','#e31a1c',
                        '#ff7f00','#cab2d6'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,10000))

axis(2,at=seq(0,10000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=5000,cex=2)
mtext('Email',1,line=1.5,at=coords[1],cex=2)
mtext('Online portal',1,line=1.5,at=coords[2],cex=2)
mtext('Phone',1,line=1.5,at=coords[3],cex=2)
mtext('Text',1,line=1.5,at=coords[4],cex=2)
mtext('Other',1,line=1.5,at=coords[5],cex=2)
dev.off()


# Bar graph 18: Share with MD
pdf(file='shared_with_md.pdf',height=6,width=5,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(1539,298,4354),
                  col=c('#4292c6','darkgray','#ef6548'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,5000))

axis(2,at=seq(0,5000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=2500,cex=2)
mtext('Yes',1,line=1.5,at=coords[1],cex=2)
mtext('Unsure',1,line=1.5,at=coords[2],cex=2)
mtext('No',1,line=1.5,at=coords[3],cex=2)
dev.off()

# Bar graph 19: Why not use
pdf(file='why_not.pdf',height=5,width=8,
    pointsize=5)
par(oma=c(1,1,1,1),mar=c(7,7,5,1))

coords <- barplot(c(41,9398,5040,584,150),
                  col=c('#1f78b4','#33a02c','#e31a1c',
                        '#ff7f00','#cab2d6'),
                  yaxt='n',xaxt='n',bty='n',border = NA,ylim=c(0,10000))

axis(2,at=seq(0,10000,1000),las=2,cex.axis=1.8)
mtext('Count',2,line=6,at=5000,cex=2)
mtext('Privacy',1,line=1.5,at=coords[1],cex=2)
mtext('Online portal',1,line=1.5,at=coords[2],cex=2)
mtext('Phone',1,line=1.5,at=coords[3],cex=2)
mtext('Text',1,line=1.5,at=coords[4],cex=2)
mtext('Other',1,line=1.5,at=coords[5],cex=2)
dev.off()