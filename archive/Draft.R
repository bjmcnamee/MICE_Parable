library("mice")
library("readxl")
library(ggplot2)
library(reshape2)
set.seed(500)

##### CLEANING #####
working_folder <- '~/rstudio/scripts/Parable/' # edit this folder to match your working folder
source(paste(working_folder,'functions.R',sep=''))
data <- data.frame(read_xlsx(paste(working_folder,'/input/PARABLE_FinalForMICE_Study_Database_31.5.2021.xlsx',sep='')))

data <- subset_LAVI_and_BNP_data(data)
data <- clean_data(data)
data <- drop_implied_features(data)
data <- rename_features(data)
features_NA <- names(data[,sapply(data, function(x) sum(is.na(x))>0)])
mice_rejected <- c('OTHER_DIURETIC_BL', 'X9M_ABPM_SYSTOLIC_24H')
features_NA <- sort(features_NA[!features_NA %in% mice_rejected])


# Run MICE on full dataset
# imp_qp1 <- mice(data, pred=quickpred(data, mincor=.3), m=10, maxit=10, meth='pmm') # quick prediction mode with minimum correlation 30% if required
# log_qp1 <- imp_qp1$loggedEvents
imp1 <- mice(data, m=10, maxit=10, meth='pmm')
log1 <- imp1$loggedEvents
completedData <- list()
for (i in 1:10) {
  completedData[i] <- list(complete(imp1,i))
  write.csv(completedData[i], (paste(working_folder,'completedData',i,'.csv',sep='')), row.names=FALSE)
}


##### REGRESSIONS #####
quit <- FALSE; pop_feat <- FALSE; loop <- 0; one_more=FALSE
while (quit==FALSE) {
  loop <- loop + 1; f_count<- 1
  ac_all1 <- sum(is.na(data))
  for (feature in (features_NA)){
    print(feature)
    fc_before <- sum(is.na(data[feature]))
    if (fc_before>0) {
      ac_before1 <- sum(is.na(data))
      print(paste('Outer:', feature, ifelse(pop_feat==TRUE,paste(':: POP feature Loop'),''), ':: LOOP',loop,'Fwd :: NA count #',fc_before,'/',ac_before1))
      test_before1 <- as.numeric(rownames(data[is.na(data[feature]),]))
      # get test features ranked and coefficients
      test_features_df <- get_features(reg_data[[f_count]], data)
      test_features <- unname(c(unlist(test_features_df[1])))
      coefficients <- unname(c(unlist(test_features_df[2])))
      data <- apply_regression(test_features, coefficients, feature, data, pop_feat)
      
      test_after <- as.numeric(rownames(data[is.na(data[feature]),]))
      ac_after1 <- sum(is.na(data))
      print(paste('Outer: End Total NA count #',ac_after1))}
    else
    {print(paste('Outer:', feature,'No NAs found - SKIPPING'))}
    print('----------------')
    f_count<- f_count+ 1
  }
  ac_all2 <- sum(is.na(data))
  features_NA = rev(features_NA); print('REVERSING FEATURES ORDER...')
  if (ac_all1==ac_all2 & one_more==TRUE) {pop_feat=TRUE; print('NONE FOUND, TRYING POP FEATURE...')}
  if (ac_all1==ac_all2) {one_more==TRUE; print('NONE FOUND, TRYING P1 MORE LOOP...')}
  if (loop==3) {quit=TRUE; print('Loop = 3 QUIT')}
}


# create regression variable list
features_NA <- names(data[,sapply(data, function(x) sum(is.na(x))>0)])
mice_rejected <- c('OTHER_DIURETIC_BL', 'X9M_ABPM_SYSTOLIC_24H')
features_NA <- features_NA[!features_NA %in% mice_rejected]
reg_data <- list()
# pool coefficients and standard errors across all regression models
# In reg_data2[i] <- find_regressors(imp2, x, y) : number of items to replace is not a multiple of replacement length
for (feature in features_NA) {
  x <- features_NA[!features_NA %in% feature]
  reg_data[[length(reg_data) + 1]] <- find_regressors(imp2, x, feature)
}
names(reg_data) <- features_NA
# reg_databkp <- reg_data



print(paste(sum(is.na(data_lean[,endpoints[1]])), sum(is.na(data_lean[,endpoints[2]])))) # count NAs in endpoint columns before 
data_lean[as.numeric(rownames(df_end_fake)),endpoints]<-NA # change same selection to NAs in main dataset
print(paste(sum(is.na(data_lean[,endpoints[1]])), sum(is.na(data_lean[,endpoints[2]])))) # count NAs in endpoint columns after

imp_qp5 <- mice(data_lean, pred=quickpred(data_lean, mincor=.3), m=10, maxit=10, meth='pmm')
imp5 <- mice(data_lean, m=10, maxit=10, meth='pmm')
log_qp5 <- imp_qp5$loggedEvents
log5 <- imp5$loggedEvents
completedData5 <- complete(imp5,1)

# compare with imputed values with true values
completedData
all <- completedData[as.numeric(rownames(df_end_fake)),endpoints]
names(all)[1] <- 'Lavi_Imputed'
names(all)[2] <- 'NTPROBNP_Imputed'
all <- cbind(all,df_end_true)
names(all)[3] <- 'Lavi_True'
names(all)[4] <- 'NTPROBNP_True'
all_NTPROBNP <- cbind(all[2],all[4])
# diff_NTPROBNP <- round((all_NTPROBNP[2]-all_NTPROBNP[1])/all_NTPROBNP[2]*100,0) # calculate % differences between imputed and true
all_LAVI <- cbind(all[1],all[3])
all_NTPROBNP <- melt(as.matrix(all_NTPROBNP))
all_LAVI <- melt(as.matrix(all_LAVI))
ggplot(all_LAVI, aes(Var1, value)) +  geom_line(aes(colour = Var2)) + ylim(20,85) # range(df_end,na.rm=T)
ggplot(all_NTPROBNP, aes(Var1, value)) +  geom_line(aes(colour = Var2)) + ylim(16,2367)




summary(imp1)
log1
summary(imp_qp1)
log_qp1
summary(imp2)
log2
summary(imp_qp2)
log_qp2
summary(imp3)
unique(log3$out)
summary(imp_qp3)
log_qp3
summary(imp4)
log4
summary(imp_qp4)
log_qp4

write.csv(data, (paste(working_folder,'output/PARABLE_cleaned.csv',sep='')), row.names=FALSE)
##### CLEANING #####




# assign variables to groups
all_cols <- names(data)
endpoints_LAVI_MAX <- get_endpoints_LAVI(all_cols)
endpoints_NTPROBNP <- get_endpoints_NTPROBNP(all_cols)
final_endpoints <- get_final_endpoints(all_cols)
measures_LAVI_MAX <- get_measures_LAVI(all_cols)
measures_NTPROBNP <- get_measures_NTPROBNP(all_cols)
endpoints_and_measures <- get_endpoints_and_measures(all_cols)
covariates <- get_covariates(all_cols)


col_index <- grep(feature, names(temp_df)) # find column index/nr of feature
temp_df <- cbind(temp_df[-col_index], temp_df[col_index]) # swap last column with feature column


# split endpoints by measurement type/name
df_end <- data.frame()
# concern wrt extra/missing variables in Sc
# placeholders - DaysSinceLastVisitSc may be 0 but other endpoints were not measured
endpoints_Sc <- c(c('X9SFEnergySc','X9SFFulloflifeSc','X9SFTiredSc','X9SFWornOutSc'),endpoints_Sc,c('DaysSinceLastVisitSc'))
# Hgt does not appear in others, assume this is not really an endpoint measurement, maybe height?
endpoints_Sc <- endpoints_Sc[!endpoints_Sc %in% c('HgtSc')] 
df_end <- rbind(df_end, data.frame(endpoints_Sc))
df_end <- cbind(df_end, data.frame(endpoints_V2))
df_end <- cbind(df_end, data.frame(endpoints_V3))
df_end <- cbind(df_end, data.frame(endpoints_V4))
df_end

# COMPLETE CASES ONLY
data_cc <- data_test[complete.cases(data_test),]
dim(data_cc)
endpoints_change <- data_cc[endpoints_V4] - data_cc[endpoints_V2]
names(endpoints_change) <- gsub('V4','_V4V2', names(endpoints_change))
data_cc <- cbind(data_cc, endpoints_change) # add to record dataframe
vars <- gsub('V4','',endpoints_V4)
model <- matrix(ncol = 0, nrow = 16)
for (var in vars) {
  col <- grep(paste(var,'_V4V2',sep=''), names(data_cc)) # get col index of 1st change variable
  fit_var <- list()
  # lm (y ~ x) - lm(Ferritin ~ DoseGroupNo + Age + SmokingCurrent + AlcoholUnitsPerWeek + GSRSOtherProds + HMPLength + HMPsoakthrough + HMPdoubleprotection + HMPLargeClots + HMPInterfereReglifestyle + HMPPeriodPain + HMPEasyBruising + HMPMedications + GSRSOtherProdsCheck)
  x <- data_cc[c(treatment, covariates)]
  y <- unlist(data_cc[col])
  fit_var <- lm(y ~ ., data = x)
  model <- cbind(model, data.frame(fit_var[1]))
  names(model)[ncol(model)] <- var
}
model
Treatment <- data_cc[,'DoseGroupNo']; Ferritin <- data_cc[,'Ferritin_V4V2']; Iron <- data_cc[,'Iron_V4V2']
ggplot(plot_df, aes(Treatment, Ferritin, fill=factor(Dose))) + geom_boxplot() + guides(fill=guide_legend(title='Treatment Dose (mg)'))
ggplot(plot_df, aes(Treatment, Iron, fill=factor(Dose))) + geom_boxplot() + guides(fill=guide_legend(title='Treatment Dose (mg)'))

# correlation
data <- data_cc
is_fac <- c(names(data[,sapply(data, is.factor)]))
df_num <- data_cc[!names(data) %in% is_fac]
#df_cor <- data.frame(cor(df_num, method = "pearson"))
#dim(df_cor)
corr_cross(df_num, # name of dataset
           max_pvalue = 0.05, # display only significant correlations (at 5% level)
           top = 25 # display top 10 couples of variables (by correlation coefficient)
)
Ferr <- df_num['Ferritin_V4V2']
corr_var(df_num, # name of dataset
         Ferr, # name of variable to focus on
         top = 25 # display top 5 correlations
) 

# T-TEST
df_t_test <- data.frame(nrow=0, ncol=0) # creates placeholder row
for (var in vars) {
  before <- data_cc[,paste(var,'V2',sep='')]
  after <- data_cc[,paste(var,'V4',sep='')]
  test <- t.test(before, after, paired=TRUE)
  df_t_test <- rbind(df_t_test, c(var,test$p.value))
}
names(df_t_test)[1] <- 'variable';names(df_t_test)[2] <- 'p-value'
df_t_test <- df_t_test[-c(1),] # delete placeholder row
df_t_test <- df_t_test[order(df_t_test[,2]),]
df_t_test
#format(test$p.value, scientific = F)
fac <- c(names(data_test[,sapply(data_cc, is.factor)]))
data_cc_no_fac <- data_cc[!names(data_cc) %in% fac]
ttest <- data.frame()


# COMPLETE + IMPUTATIONS
# create imputations
inlist <- c(covariates, treatment)
data_test_in <- data_test
# quickpred
impq <- mice(data_test_in, pred=quickpred(data_test, mincor=.3), seed=500, m=1, maxit=1, meth='pmm', inlist=inlist) # save imputations
logq <- impq$loggedEvents
outlist2 <- unique(gsub(' ','',unlist(strsplit(paste(logq$out, collapse=','),','))))
outlist2
length(outlist2)
dim(data_test)
data_test_in <- data_test[!names(data_test) %in% outlist2]
dim(data_test_in)
impq <- mice(data_test_in, pred=quickpred(data_test_in, mincor=.3), seed=500, m=1, maxit=1, meth='pmm', inlist=inlist) # save imputations
logq <- impq$loggedEvents
logq
# all/standard
impq <- mice(data_test_in, seed=500, m=1, maxit=1, meth='pmm', inlist=inlist) # save imputations
logq <- impq$loggedEvents
outlist <- unique(gsub(' ','',unlist(strsplit(paste(logq$out, collapse=','),','))))
outlist
length(outlist)
dim(data_test)

pm <- names(data.frame(imp$predictorMatrix))
setdiff(pm, covariates)
setdiff(covariates, pm)
covariates
pm

completedData <- list() # create empty list
data_test_imp <- data_test_in
# loop through imp - save imputations + true values to data_test
for (i in 1:m) {
  length(completedData[i])
  completedData[i] <- list(complete(imp,i)) # add dataset to completedData list
  imputations <- data.frame(completedData[i])[endpoints_V4] - data.frame(completedData[i])[endpoints_V2]
  names(imputations) <- gsub('V4','', names(imputations))
  names(imputations) <- c(paste(names(imputations),'V4V2',i,sep='_'))
  data_test_imp <- cbind(data_test_imp, imputations) # add to record dataframe
}
dim(data_test_imp)
names(data_test_imp)


# analyse pm before removing variables
pm <- data.frame(rowSums(impq$predictorMatrix))
names(pm)[1] <- 'predictor'
pm_not <- subset(pm, pm[,1]==0)
pm_not <- c(rownames(pm_not))
pm_not
pm <- subset(pm, pm[,1]>0)
pm <- rownames(pm)
pm


pred_matrix <- make.predictorMatrix(data_test_imp, blocks = make.blocks(data_test_imp))
ncol(pred_matrix)
View(pred_matrix)
pred_matrix <- pred_matrix[c(treatment, covariates, endpoints)]

# all
df_end
i <- 9
endpoint <- c(unlist(strsplit(apply(df_end, 1, paste, collapse=',')[i],',')))
all_vars <- c(treatment, covariates, endpoint)
all_vars
names(data_test[all_vars])
length(names(data_test[all_vars]))
data <- data_test[all_vars]

m <- 2 # datasets
imp <- mice(data, seed=500, m=m, maxit=1, meth='pmm') # save imputations
log <- imp$loggedEvents
log
completedData <- list() # create empty list
data_test_imp <- data_test
# loop through imp - save imputations + true values to data_test
for (i in 1:m) {
  completedData[i] <- list(complete(imp,i)) # add dataset to completedData list
  imputations <- data.frame(completedData[i])[endpoints_V4] - data.frame(completedData[i])[endpoints_V2]
  names(imputations) <- gsub('V4','', names(imputations))
  names(imputations) <- c(paste(names(imputations),'V4V2',i,sep='_'))
  data_test_imp <- cbind(data_test_imp, imputations) # add to record dataframe
}
dim(data_test_imp)
names(data_test_imp)

# loop through all imp cols, run lm for each 
var <- 'Ferritin'
var <- 'Iron'
col <- grep(paste(var,'_V4V2_1',sep=''), names(data_test_imp)) # get col index of 1st change variable
fit_var <- list()
for (i in 1:m) {
  # lm (y ~ x) - lm(Ferritin ~ DoseGroupNo + Age + SmokingCurrent + AlcoholUnitsPerWeek + GSRSOtherProds + HMPLength + HMPsoakthrough + HMPdoubleprotection + HMPLargeClots + HMPInterfereReglifestyle + HMPPeriodPain + HMPEasyBruising + HMPMedications + GSRSOtherProdsCheck)
  x <- data_test_imp[c(treatment, covariates)]
  y <- unlist(data_test_imp[col])
  fit_var[i] <- list(lm(y ~ ., data = x))
  col <- col + 31
}
pool(fit_var)
summary <- summary(pool(fit_var))
summary






############# END #############
# find factors
is_fac <- c(names(data_test[,sapply(data_test, is.factor)]))
is_fac
df_num <- data_test[!names(data_test) %in% is_fac]
# find cor
cor_relation <- cor(df_num)
cor_relation[abs(cor_relation) < 0.7] <- NA
cor_relation

names(data_test[sapply(data_test, function(x) length(unique(x))) == nrow(data_test)-1]) # which cols have only unique values
# keep DOB and remove other columns
DOB <- data_test['DOB']
data_test <- Filter(function(x) length(unique(x)) != nrow(data_test)-1, data_test)
data_test <- cbind(data_test, DOB)

# contains study notes and section descriptions (section is several related columns) with no data - imported as column names
names(Filter(function(x)all(is.na(x)), data_test)) # 26 count
# remove empty columns
data_test <- Filter(function(x)!all(is.na(x)), data_test)
dim(data_test) # rows x cols

dim(data_test) # rows x cols
names(data_test) # names of variables/cols
head(data_test) # peak at data - top five rows
str(data_test, list.len=ncol(data_test)) # data stucture - object (eg dataframe), dimensions, variable types (eg factor, numeric, character, etc)

# lm standard
dx <- data_test[x]
idx <- sapply(dx, is.factor)
dx[idx] <- lapply(dx[idx], function(x) as.numeric(as.character(x)))
str(dx)
x <- dx[idx]
fit_standard <- with(imp, lm(as.formula(paste('scale(',y,')~scale(',paste(x,collapse =' + '),')'))))

# rename in dataframe analysis_vars
names(data_test)[names(data_test) == 'DaysSinceSc'] <- 'DaysSinceLastVisitV2'
names(data_test)[names(data_test) == 'DaysSinceV2'] <- 'DaysSinceLastVisitV3'
names(data_test)[names(data_test) == 'DaysSinceV3'] <- 'DaysSinceLastVisitV4'

######## END ############

# previously used code replaced by two van Burren lines above

completedData <- list() # create empty list
# loop through imp - save imputations + true values (only 'JDMRI_LAVI_MAX_BL','CHGLAVImax') to data_test
for (i in 1:50) {
  completedData[i] <- list(complete(imp,i)) # add dataset to completedData list
  imputations <- data.frame(completedData[i])[c('JDMRI_LAVI_MAX_BL','CHGLAVImax')]
  names(imputations) <- c(paste('imp_BL',i,sep='_'),paste('imp_CHG',i,sep='_'))
  data_test <- cbind(data_test, imputations) # add to record dataframe
}
# loop through 50 cols ('imp_CHG_1'...'imp_CHG_50'), run lm for each 
col <- 10 # counter for 'imp_CHG_X' col
fitLAVImax <- list()
for (i in 1:50) {
  # lm (y ~ x) - lm(CHGLAVImax ~ TREATMENT_ALLOCATION + AGE_BASELINE + GENDER + VASC_BL + HPT_BL + DIABETES_BL + OBESE_BL + JDMRI_LAVI_MAX_BL)
  x <- data_test[c(1:7,col)] # JDMRI_LAVI_MAX_BL = 'imp_BL_X' for x 1 to 10
  y <- unlist(data_test[col+1]) # CHGLAVImax = 'imp_CHG_X' for x 1 to 10
  fitLAVImax[i] <- list(lm(y ~ ., data = x))
  col <- col + 2
}
pool(fitLAVImax)
summary <- summary(pool(fitLAVImax))[1:8,]


library('stringr')
library('dplyr')
library('zoo')
#### FUNCTIONS #####

covariates <- c('TREATMENT_ALLOCATION', 'AGE_BASELINE', 'GENDER', 'BL_TG', 'BL_LDLC', 'BL_SBP', 'BL_DBP', 'BL_PULSE', 'MACECOUNT_POST', 'CAD_BL', 'ANGINA_BL', 'MI_BL', 'IHD_BL', 'HPT_BL', 'DYSLIPIDEMIA_BL', 'DIABETES_BL', 'STROKE_BL', 'TIA_BL', 'ANTIARRTHYMIC_BL', 'ALPHABLOCKER_BL', 'BETABLOCKER_BL', 'CCB_BL', 'ALDOSTERONE_ANTAG_BL', 'STATIN_BL', 'OTHER_DYSLIPID_BL', 'THIAZIDE_BL', 'LOOP_BL', 'OTHER_DIURETIC_BL', 'ANTIPLATELET_EXCL_ASP_BL', 'ASPIRIN_BL', 'ANTICOAG_EXCL_WARFARIN_BL', 'WARFARIN_BL', 'INSULIN_BL', 'METANTIDIABETIC_BL', 'SUREAANTIDIABETIC_BL', 'DPP4ANTIDIABETIC_BL', 'GLPANTIDIABETIC_BL', 'SGLT2ANTIDIABETIC_BL', 'BMI_BL')
baseline <- c('BL_NTPROBNP','BL_BNP', 'BL_ABPM_SYSTOLIC_24H', 'BL_ABPM_DIASTOLIC_24H', 'BL_ABPM_PULSEPRESSURE_24H', 'BL_ABPM_HEARTRATE_24H', 'ECHO_LAVI_BL', 'JDMRI_LAVI_MAX_BL', 'ECHO_E_Ep_AVG_BL')

subset_LAVI_and_BNP_data <- function(data) {
  print(paste('Dataset dimensions (rows x cols):', nrow(data),'x ',ncol(data)))
  print('Reducing full dataset to covariates and NTPROBNP & LAVI endpoints and measurements...')
  all_variables <- c('TREATMENT_ALLOCATION', 'AGE_BASELINE', 'GENDER', 'BL_TG', 'BL_LDLC', 'MACECOUNT_POST', 'CAD_BL', 'ANGINA_BL', 'MI_BL', 'IHD_BL', 'HPT_BL', 'DYSLIPIDEMIA_BL', 'DIABETES_BL', 'STROKE_BL', 'TIA_BL', 'ANTIARRTHYMIC_BL', 'ALPHABLOCKER_BL', 'BETABLOCKER_BL', 'CCB_BL', 'ALDOSTERONE_ANTAG_BL', 'STATIN_BL', 'OTHER_DYSLIPID_BL', 'THIAZIDE_BL', 'LOOP_BL', 'OTHER_DIURETIC_BL', 'ANTIPLATELET_EXCL_ASP_BL', 'ASPIRIN_BL', 'ANTICOAG_EXCL_WARFARIN_BL', 'WARFARIN_BL', 'INSULIN_BL', 'METANTIDIABETIC_BL', 'SUREAANTIDIABETIC_BL', 'DPP4ANTIDIABETIC_BL', 'GLPANTIDIABETIC_BL', 'SGLT2ANTIDIABETIC_BL', 'BMI_BL', 'BL_SBP', 'BL_DBP', 'BL_PULSE', 'JDMRI_LAVI_MAX_BL', 'JDMRI_LAVImax_18M', 'BL_NTPROBNP', 'X3M_NTPROBNP', 'X6M_NTPROBNP', 'X9M_NTPROBNP', 'X12M_NTPROBNP', 'X15M_NTPROBNP', 'X18M_NTPROBNP', 'BL_BNP', 'X3M_BNP', 'X6M_BNP', 'X9M_BNP', 'X12M_BNP', 'X15M_BNP', 'X18M_BNP', 'BNPchg9M', 'BNPchg18M', 'BL_ABPM_SYSTOLIC_24H', 'BL_ABPM_DIASTOLIC_24H', 'BL_ABPM_PULSEPRESSURE_24H', 'BL_ABPM_HEARTRATE_24H', 'X9M_ABPM_SYSTOLIC_24H', 'X9M_ABPM_DIASTOLIC_24H', 'X9M_ABPM_PULSEPRESSURE_24H', 'X9M_ABPM_HEARTRATE_24H', 'X18M_ABPM_SYSTOLIC_24H', 'X18M_ABPM_DIASTOLIC_24H', 'X18M_ABPM_PULSEPRESSURE_24H', 'X18M_ABPM_HEARTRATE_24H', 'ECHO_LAVI_BL', 'ECHO_E_LAT_Ep_BL', 'ECHO_LAVI_9M', 'ECHO_E_Ep_AVG_9M', 'ECHO_LAVI_18M', 'ECHO_E_Ep_AVG_18M', 'CHGLVEF18M', 'CHGLVMI18M', 'ChangeJDMRI_LA_EF_18M', 'CHGLASVI')
  data <- data[,names(data) %in% all_variables] 
  print(paste('Dataset dimensions (rows x cols):', nrow(data),'x',ncol(data)))
  return(data)
}


drop_excess_NA_features <- function(data, col_threshold, row_threshold, keep_endpoints) {
  print('Removing features with excess NAs...')
  endpoints_LAVI <- names(data[,grep('JDMRI_LAVI_MAX',names(data))])
  endpoints_NTPROBNP <- names(data[,grep('NTPROBNP',names(data))])
  endpoints <- c(endpoints_LAVI, endpoints_NTPROBNP)
  # remove cols with > X% missing 
  col_missing <- unlist(sapply(data, function(x) (sum(is.na(x))/length(x)*100)))
  all_cols <- ncol(data)
  drop_cols <- names(data[col_missing>col_threshold])
  drop_cols_some <- drop_cols[!drop_cols %in% endpoints]
  if (keep_endpoints==TRUE) {
    data <- data[!names(data) %in% drop_cols_some]
    sub_cols <- ncol(data)} 
  else {
    data <- data[!names(data) %in% drop_cols]
    sub_cols <- ncol(data)}
  cols_removed <- all_cols - sub_cols
  print(paste('Dropped', cols_removed,'columns with >',col_threshold,'% missing values :'))
  print(cols_removed)
  # remove rows with > X% missing 
  row_missing <- apply(data, 1, function(x) round((sum(is.na(x))/(sum(is.na(x))+sum(!is.na(x))))*100,0))
  all_rows <- nrow(data)
  rows_removed <- as.numeric(rownames(data[row_missing > row_threshold,]))
  data <- data[row_missing < row_threshold,]
  sub_rows <- nrow(data)
  print(paste('Dropped', all_rows-sub_rows,'observations with >',row_threshold,'% missing values :'))
  print(rows_removed)
  print(paste('Dataset dimensions (rows x cols):', nrow(data),'x',ncol(data)))
  return(data)
}

drop_implied_features <- function(data) {
  print('Removing features with implied values* :')
  implied_features <- c('BNPchg9M', 'BNPchg18M', 'CHGLASVI', 'CHGLVEF18M', 'CHGLVMI18M', 'ChangeJDMRI_LA_EF_18M')
  print(implied_features)
  print('*Tested BNPchg9M = X9M_BNP - BL_BNP (note 13/250 discrepancies)')
  data <- data[,!names(data) %in% implied_features]
  print(paste('Dataset dimensions (rows x cols):', nrow(data),'x',ncol(data)))
  return(data)
}

rename_features <- function(data) {
  # rename columns with'.' to'_'
  print('Renaming features with . in column name : . to _')
  names(data) <- gsub(x = names(data), pattern = '\\.', replacement = '_')
  # rename features with time suffix - standardise naming convention
  print('Renaming features with time suffix - standardise naming convention')
  names(data)[names(data) == 'ECHO_LAVI_BL'] <- 'BL_ECHO_LAVI'
  names(data)[names(data) == 'ECHO_E_Ep_AVG_BL'] <- 'BL_ECHO_E_Ep_AVG'
  names(data)[names(data) == 'ECHO_LAVI_9M'] <- 'X9M_ECHO_LAVI'
  names(data)[names(data) == 'ECHO_E_Ep_AVG_9M'] <- 'X9M_ECHO_E_Ep_AVG'
  names(data)[names(data) == 'ECHO_LAVI_18M'] <- 'X18M_ECHO_LAVI'
  names(data)[names(data) == 'ECHO_E_Ep_AVG_18M'] <- 'X18M_ECHO_E_Ep_AVG'
  names(data)[names(data) == 'JDMRI_LAVI_MAX_BL'] <- 'BL_JDMRI_LAVI_MAX'
  names(data)[names(data) == 'JDMRI_LAVImax_18M'] <- 'X18M_JDMRI_LAVI_MAX' # JDMRI_LAVImax_18M has lowercase aswell as suffix
  print(paste('Dataset dimensions (rows x cols):', nrow(data),'x',ncol(data)))
  return(data)
}

get_endpoints_LAVI <- function(cols) {
  endpoints_LAVI <- sort(c(cols[grep('JDMRI_LAVI_MAX',cols)]))
  #print(paste('Endpoints LAVI_MAX x',length(endpoints_LAVI),':',toString(endpoints_LAVI)))
  return(endpoints_LAVI)
}
get_endpoints_NTPROBNP <- function(cols) {
  endpoints_NTPROBNP <- sort(c(cols[grep('NTPROBNP',cols)]))
  #print(paste('Endpoints NTPROBNP x',length(endpoints_NTPROBNP),':',toString(endpoints_NTPROBNP)))
  return(endpoints_NTPROBNP)
}
get_final_endpoints <- function(cols) {
  final_endpoints <- sort(c(endpoints[grep('18M',endpoints)]))
  #print(paste('Endpoints - Final x',length(final_endpoints),':',toString(final_endpoints)))
  return(final_endpoints)
}
get_measures_LAVI <- function(cols) {
  measures_LAVI <- sort(c(cols[grep('ECHO',cols)]))
  #print(paste('Measures LAVI_MAX x',length(measures_LAVI),':',toString(measures_LAVI)))
  return(measures_LAVI)
}
get_measures_NTPROBNP <- function(cols) {
  measures_NTPROBNP <- sort(c(cols[grep('_BNP|_ABPM',cols)]))
  #print(paste('Measures BNP, ABPM x',length(measures_NTPROBNP),':',toString(measures_NTPROBNP)))
  return(measures_NTPROBNP)
}
get_endpoints_and_measures <- function(cols) {
  endpoints_and_measures <- sort(c(get_endpoints_LAVI(cols), get_endpoints_NTPROBNP(cols), get_measures_LAVI(cols), get_measures_NTPROBNP(cols)))
  #print(paste('Endpoints and Measures x',length(endpoints_and_measures),':',toString(endpoints_and_measures)))
  return(endpoints_and_measures)
}
get_covariates <- function(cols) {
  covariates <- sort(c(cols[!cols %in% endpoints_and_measures]))
  #print(paste('Covariates x',length(covariates),':',toString(covariates)))
  return(covariates)
}

IsOutlier <- function(data_tmp) {
  lowerq = quantile(data_tmp, na.rm = TRUE)[2]
  upperq = quantile(data_tmp, na.rm = TRUE)[4]
  iqr = upperq - lowerq 
  threshold_upper = (iqr * 1.5 * 2) + upperq
  threshold_lower = lowerq - (iqr * 1.5 * 2)
  data_tmp > threshold_upper | data_tmp < threshold_lower 
}

remove_outliers <- function(data_tmp, threshold) {
  data_factors <- cbind(data_tmp[1],data_tmp[3])
  data_tmp <- cbind(data_tmp[2],data_tmp[4:length(data_tmp)])
  row_outliers <- as.numeric(rownames(data_tmp[rowSums(sapply(data_tmp, IsOutlier), na.rm = TRUE) > threshold,]))
  data_tmp <- cbind(data_factors,data_tmp)
  data_tmp <- data_tmp[-row_outliers, ]
  return(data_tmp)
} 

get_t_test_data <- function(n, var1, var2) {
  test <- t.test(var1, var2, paired=TRUE)
  # variable name
  var_name <- paste(n,'_',toupper(substr(deparse(substitute(var1)),1,2)),'-',toupper(substr(deparse(substitute(var2)),1,2)),sep='')
  # t test results
  t_test_significant <- ifelse(test$p.value<0.025,'significant','')
  t_test_p_value <- round(test$p.value,4)
  t_test_test_statistic <- round(test$statistic[[1]],4)
  t_test_confidence_interval_lower <- round(test$conf.int[1],4)
  t_test_confidence_interval_upper <- round(test$conf.int[2],4)
  t_test_standard_error <- round(test$stderr,4)
  # Shapiro Wilks normality test
  diff <- var1 - var2
  shap.test <- shapiro.test(diff)
  normality <- ifelse(shap.test$p.value<0.05,'Normal','')
  shapiro_p.value <- round(shap.test$p.value,4)
  outliers = length(boxplot(diff)$out)
  # add values to vector
  test <- c(var_name, t_test_significant, t_test_p_value, t_test_test_statistic, t_test_confidence_interval_lower, t_test_confidence_interval_upper, t_test_standard_error, normality, shapiro_p.value, outliers)
  return(test)
}

get_t_test_report_active_iron <- function(data) {
  # subset V2, V3, V4 column name list s
  v2 <- data[ , grepl('V2' , names(data) ) ]
  v3 <- data[ , grepl('V3' , names(data) ) ]
  v4 <- data[ , grepl('V4' , names(data) ) ]
  # remove V2, V3, V4 from column names
  names(v2) <- sort(gsub(x = names(v2), pattern ='V2', replacement =''))
  names(v3) <- sort(gsub(x = names(v3), pattern ='V3', replacement =''))
  names(v4) <- sort(gsub(x = names(v4), pattern ='V4', replacement =''))
  # create list of V2, V3, V4 lists values
  vx <- NULL
  var_all <- NULL
  for (n in names(v2)) {
    vx[[n]] <- list(v2[,n],v3[,n],v4[,n])
    var_all <- c(var_all,vx[n]) 
  }
  # run t test for all variables and save to dataframe
  test_data <- data.frame()
  for (n in names(v2)) {
    v2_var <- unlist(var_all[n][[1]][1])
    v3_var <- unlist(var_all[n][[1]][2])
    v4_var <- unlist(var_all[n][[1]][3])
    # test differences between V2, V3
    test <- get_t_test_data(n, v2_var, v3_var)
    test_data <- rbind(test_data, as.data.frame(t(test)))
    # test differences between V3, V4
    test <- get_t_test_data(n, v3_var, v4_var)
    test_data <- rbind(test_data, as.data.frame(t(test)))
    # test differences between V2, V4
    test <- get_t_test_data(n, v2_var, v4_var)
    test_data <- rbind(test_data, as.data.frame(t(test)))
  }
  names(test_data)<-c('variable','t test','p.value','statistic','conf.int1','conf.int2','stderr','Shapiro','p.value','outliers')
  write.csv(test_data,'~/rstudio/scripts/solvitron/t.test_data.csv', row.names=FALSE)
  write.csv(test_data,'~/0python/solvitron/t.test_data.csv', row.names=FALSE)
  return (test_data)
}

locf_update_data <- function(data) {
  # Last observation carried forward (LOCF)
  # create new dataframe with measures
  print('Creating new subset dataframe with measures and endpoints only')
  all_cols <- names(data)
  endpoints <- c(get_endpoints_LAVI(all_cols),get_endpoints_NTPROBNP(all_cols))
  measures <- c(get_measures_LAVI(all_cols),get_measures_NTPROBNP(all_cols))
  data_locf <- data[, c(measures, endpoints)] # subset time repeating group columns
  # count NAs before LOCF
  before <- sum(is.na(data_locf))
  # create a copy to compare results
  write.csv(data_locf,'~/rstudio/scripts/Parable/output/data_locf.before.csv', row.names=FALSE)
  # create list of distinct variables measured
  print('Found following groups')
  groups <- c(unique(str_split(names(data_locf),'_', n=2, simplify=T)[,2]))
  print(toString(groups))
  # remove groups with only 1 element - na.locf function will fail
  for (group in groups[1:length(groups)]) {
    if (length(grep(group, names(data_locf))) < 2) {groups = groups[!(groups %in% group)]}
  }
  # loop through each group in groups list, apply na.locf function to each column in data_locf
  data_locf.new <- data.frame(matrix(0, ncol = 0, nrow = 250)) # create empty df
  print(paste('Processing groups...'))
  for (group in groups[1:length(groups)]) {
    # subset the data_locf frame by indexing on the cols_locf columns (3rd column) matching the group and save to temp dataframe
    group_cols <- grep(group, names(data_locf))
    print(paste(group,':', toString(names(data_locf[,group_cols]))))
    data_locf.tmp <- data_locf[, group_cols]
    # apply the na.locf function to temp dataframe
    data_locf.tmp <- data.frame(t(apply(data_locf.tmp, 1, function(x) na.locf(x, fromLast = F, na.rm = F))))
    data_locf.new <- cbind(data_locf.new, data_locf.tmp)
  }
  # count NAs after LOCF
  after <- sum(is.na(data_locf.new))
  # calculate % difference
  reduction <- round((before-after)/before*100,1)
  # display result
  print(paste('NAs:: Before:',before,'After:',after,'Reduction(%):',reduction))
  write.csv(data_locf.new,'~/rstudio/scripts/Parable/output/data.new.locf.csv', row.names=FALSE)
  return(data_locf.new)
}

find_regressors <- function(imp, x, y) {
  done=FALSE
  print(paste('Running step wise regression on',y,'...'))
  while (done==FALSE) {
    fit <- with(imp1, lm(as.formula(paste(y,'~',paste(x,collapse =' +')))))
    pool_fit <- pool(fit) # pool coefficients and standard errors across all regression models
    summary <- summary(pool(fit))
    if (summary[1,][1] =='(Intercept)'){intercept <- summary[1,]}
    max_p_line <- summary[summary$p.value == max(summary$p.value), ]
    max_p_var <- toString(max_p_line[,'term'])
    if (max_p_var=='TREATMENT_ALLOCATIONb') max_p_var='TREATMENT_ALLOCATION' # factor error handling
    if (max_p_var=='GENDERmale') max_p_var='GENDER' # factor error handling
    max_p_val <- toString(max_p_line[,'p.value'])
    if (max_p_val > 0.05) {x <- x[!x %in% max_p_var]}
    else {done=TRUE}
  }
  return(summary)
}

get_features <- function(reg_data, data) {
  test_features_df <- data.frame(reg_data[1])
  test_features_df <- data.frame(lapply(test_features_df, as.character), stringsAsFactors=FALSE) # save model features
  names(test_features_df)[1] <- 'features'
  coefficients <- reg_data[2] # save model coefficients
  test_features_df <- cbind(test_features_df, coefficients)
  cols <- c(as.character(unlist(test_features_df[1])))
  means <- colMeans(data[,cols],na.rm=TRUE)
  test_features_df <- cbind(test_features_df, means)
  names(test_features_df)[3] <- 'mean'
  rownames(test_features_df)<- c()
  test_features_df <- cbind(test_features_df, test_features_df[,2]*test_features_df[,3])
  names(test_features_df)[4] <- 'product'
  test_features_df <- cbind(test_features_df, rank(abs(test_features_df[,4])))
  names(test_features_df)[5] <- 'rank'
  test_features_df <- test_features_df[order(test_features_df[,5]),]
  return(test_features_df)
}

apply_regression <- function(test_features, coefficients, feature, data, pop_feat) {
  if (pop_feat == TRUE & length(test_features)>2) {
    print(paste('Popped',head(test_features,1)))
    test_features <- tail(test_features,-1)
    print(head(test_features,-1))}
  test_features <- c(test_features, feature)
  print(test_features)
  temp_df <- data[test_features] # subset master dataset by features
  if ('TREATMENT_ALLOCATION' %in% names(temp_df)) { # convert factor to numeric if exists
    temp_df[,'TREATMENT_ALLOCATION'] <- as.numeric(temp_df[,'TREATMENT_ALLOCATION'])} 
  if ('GENDER' %in% names(temp_df)) { # convert factor to numeric if exists
    temp_df[,'GENDER'] <- as.numeric(temp_df[,'GENDER'])} 
  before <- sum(is.na(temp_df[,ncol(temp_df)])) # count NAs before imputation value replacement
  
  # subset rows with NA in column of interest & no NA in all others
  NA_found <- is.na(temp_df[,ncol(temp_df)]) & !is.na(rowSums(temp_df[, -ncol(temp_df)]))
  
  
  if (dim(temp_df[which(NA_found),])[1]!=0) {
    # temp_df[which(NA_found),][ncol(temp_df)] <- coefficients[1] # save 1st coefficient (intercept) value
    # loop through all features (other features, not last feature)
    for (i in 1:(length(test_features)-1)){
      # save product (feature coefficient x feature value) to test df subset
      x <- coefficients[i]*temp_df[which(NA_found),][,i]
      if (any(NA_found) == TRUE) {temp_df[which(NA_found),][ncol(temp_df)] <- temp_df[which(NA_found),][ncol(temp_df)] + x
      }}}
  
  
  # count NAs after imputated value replacement
  after <- sum(is.na(temp_df[,ncol(temp_df)])) 
  
  if (dim(temp_df[which(NA_found),])[1]!=0) {
    print(paste('Inner :', feature,'updated - replaced',before-after,'values'))
    data[,feature] <- temp_df[,feature] } # update NAs with imputations by replacing feature column
  else {print(paste('Inner :', feature,'not updated - No NAs found'))}
  
  return(data)
}

load_clean_data <- function(data){
  data <- subset_LAVI_and_BNP_data(data)
  data <- clean_data(data)
  data <- drop_implied_features(data)
  data <- rename_features(data)
  features_NA <- names(data[,sapply(data, function(x) sum(is.na(x))>0)])
  mice_rejected <- c('OTHER_DIURETIC_BL', 'X9M_ABPM_SYSTOLIC_24H')
  features_NA <- sort(features_NA[!features_NA %in% mice_rejected])
  return(data)
}

