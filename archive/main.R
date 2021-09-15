library("mice")
library("readxl")
working_folder <- '~/rstudio/scripts/Parable/' # edit this folder to match your working folder
source(paste(working_folder,'functions.R',sep=''))
data <- data.frame(read_xlsx(paste(working_folder,'/input/PARABLE_FinalForMICE_Study_Database_31.5.2021.xlsx',sep='')))
data_test <- data
features <- c('TREATMENT_ALLOCATION','AGE_BASELINE','GENDER','DIABETES_BL','HPT_BL','OBESE_BL','VASC_BL','JDMRI_LAVI_MAX_BL','JDMRI_LAVImax_18M')

# subset required features
data_test <- data_test[,names(data_test) %in% features] # subset data on vars
head(data_test) # peak at data - top five rows
dim(data_test) # rows x cols
str(data_test) # data stucture - object (eg dataframe), dimensions, variable types (eg factor, numeric, character, etc)
names(data_test) # names of variables/cols

# transform data
data_test <- clean_data(data_test) # basic text cleaning
test_feature <- data_test['JDMRI_LAVImax_18M']-data_test['JDMRI_LAVI_MAX_BL']
names(test_feature) <- 'CHGLAVImax'
data_test <- cbind(data_test, test_feature) # add 'CHGLAVImax' col

# create imputations
imp <- mice(data_test, seed=500, m=50, maxit=10, meth='pmm') # save imputations to imp6

# pool coefficients and standard errors across all regression models - van Burren code
x <- head(features,-1) # select x variables, all but last
y <- names(test_feature) # select y variable, ie test feature
fit <- with(imp, lm(as.formula(paste(y,'~',paste(x,collapse =' +'))))) # fit regression
pool <- pool(fit) # pool results
summary <- summary(pool) # get summary of pooled results
conf <- summary(pool, conf.int = TRUE, conf.level=0.95) # generate CIs
print(conf) # show pooled results + CIs
coef <- summary[1:2] # subset variables and coefficients
OR <- cbind(coef[1], exp(coef[2])) # add variables and ln coefficients to dataframe
OR <- cbind(OR, exp(conf[7:8])) # add ln CIs to dataframe
print(OR) # display dataframe





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
