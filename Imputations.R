library("mice")
library("readxl")
working_folder <- '~/rstudio/scripts/Parable/' # edit this folder to match your working folder
source(paste(working_folder,'functions.R',sep=''))
data <- data.frame(read_xlsx(paste(working_folder,'/input/PARABLE_FinalForMICE_Study_Database_31.5.2021.xlsx',sep='')))
treatment <- 'TREATMENT_ALLOCATION'
covariates <- c('AGE_BASELINE','GENDER','DIABETES_BL','HPT_BL','OBESE_BL','VASC_BL')
NTPROBNP <- c('BL_NTPROBNP','X3M_NTPROBNP','X6M_NTPROBNP','X9M_NTPROBNP','X12M_NTPROBNP','X15M_NTPROBNP','X18M_NTPROBNP')
JDMRI_LAVI_MAX <- c('JDMRI_LAVI_MAX_BL','JDMRI_LAVImax_18M')
ABPM_PULSEPRESSURE <- c('BL_ABPM_PULSEPRESSURE_24H','X9M_ABPM_PULSEPRESSURE_24H','X18M_ABPM_PULSEPRESSURE_24H')

# run one only of the next three lines depending on endpoint required
var_impute <- NTPROBNP
var_impute <- JDMRI_LAVI_MAX
var_impute <- ABPM_PULSEPRESSURE

name_var_impute <- gsub('^BL_|_BL$','',var_impute[1]) # get substring
features <- c(treatment, covariates, var_impute)
baseline_measure <- head(var_impute,1)
final_measure <- tail(var_impute,1)

# subset required features
data_test <- data
data_test <- data_test[,names(data_test) %in% features] # subset data on vars
data_test <- clean_data(data_test) # basic text cleaning
head(data_test) # peak at data - top five rows
dim(data_test) # rows x cols
str(data_test) # data stucture - object (eg dataframe), dimensions, variable types (eg factor, numeric, character, etc)

# create imputations
imp <- mice(data_test, seed=500, m=50, maxit=10, meth='pmm') # save imputations to imp6
imp$loggedEvents # examine log

# pool coefficients and standard errors across all regression models - van Burren code
x <- c(treatment, covariates, head(var_impute,1)) # select x variables, excluding all measures but baseline
y <- c(final_measure, baseline_measure) # select y variable, ie test feature
model <- paste(paste(y,collapse =' - '),'~',paste(x,collapse =' + '))
fit <- with(imp, lm(as.formula(model))) # fit regression
pool <- pool(fit) # pool results
summary <- summary(pool) # get summary of pooled results
conf <- summary(pool, conf.int = TRUE, conf.level=0.95) # generate CIs
conf

# Log Transform results
coef <- summary[1:2] # subset variables and coefficients
Log_Transform <- data.frame()
Log_Transform <- cbind(coef[1], sign(coef[2])*log(abs(coef[2])+1)) # add variables and ln coefficients to dataframe
Log_Transform <- cbind(Log_Transform, sign(conf[7])*log(abs(conf[7])+1)) # add ln CIs to dataframe
Log_Transform <- cbind(Log_Transform, sign(conf[8])*log(abs(conf[8])+1)) # add ln CIs to dataframe
Log_Transform

# Exp Transform results
coef <- summary[1:2] # subset variables and coefficients
Exp_Transform <- data.frame()
Exp_Transform <- cbind(coef[1], exp(coef[2])) # add variables and ln coefficients to dataframe
Exp_Transform <- cbind(Exp_Transform, exp(conf[7:8])) # add ln CIs to dataframe
Exp_Transform

write_results(name_var_impute, toString(model), conf, Log_Transform, Exp_Transform)

