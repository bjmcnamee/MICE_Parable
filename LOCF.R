library("mice")
library("readxl")
library(ggplot2)
library(reshape2)
set.seed(500)

working_folder <- '~/rstudio/scripts/Parable/' # edit this folder to match your working folder
source(paste(working_folder,'Parable_functions.R'))
data <- load_clean_data() # load raw dataset - xsheet from ML

# Last observation carried forward (LOCF) - method of imputing missing data in longitudinal studies : maintain the sample size and to reduce the bias caused by the attrition of participants in a study
test_df <- data[,which(unlist(lapply(data, function(x) any(is.na(x)))))]
names(test_df)
data_locf <- locf_update_data(data)

# plot true values versus estimated values for all available estimated values
feature <- 'X18M_JDMRI_LAVI_MAX'
plot_df <- data.frame(data[,feature],data_locf[,feature])
names(plot_df)[1] <- 'Lavi_True'
names(plot_df)[2] <- 'Lavi_Imputed'
# plot_df <- na.omit(plot_df)
plot_df <- melt(as.matrix(plot_df))
ggplot(plot_df, aes(Var1, value)) +  geom_line(aes(colour = Var2)) + ylim(20,85) # endpoint_LAVI_MAX

feature <- 'X18M_NTPROBNP'
ggplot(plot_df, aes(Var1, value)) +  geom_line(aes(colour = Var2)) + ylim(16,2367) # endpoint_NTPROBNP
