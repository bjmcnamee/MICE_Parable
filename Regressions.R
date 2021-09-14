##### REGRESSIONS #####
# create regression variable list
missing_vars <- names(data[,sapply(data, function(x) sum(is.na(x))>0)])
mice_rejected <- c('OTHER_DIURETIC_BL', 'X9M_ABPM_SYSTOLIC_24H')
missing_vars <- sort(missing_vars[!missing_vars %in% mice_rejected])
regression <- list()
# pool coefficients and standard errors across all regression models
# In regression2[i] <- find_regressors(imp2, x, y) : number of items to replace is not a multiple of replacement length
for (y in missing_vars) {
  x <- missing_vars[!missing_vars %in% y]
  regression[[length(regression) + 1]] <- find_regressors(imp2, x, y)
  }

  
source(paste(working_folder,'functions.R',sep=''))
i <- 1
for (y in missing_vars) {
  print(paste('----',y,'----'))
  if (i==1) {print(paste('NA count #',sum(is.na(data))))}
  col_index <- grep(y, names(data))
  data <- cbind(data[-col_index], data[col_index])
  data <- apply_regression(regression[[i]], data)
  print(paste('NA count #',sum(is.na(data))))
  i <- i + 1
}

# plot true values versus estimated values for all available estimated values
plot_df <- temp_df[,(ncol(temp_df)-2):(ncol(temp_df)-1)]
plot_df
head(plot_df)
dim(plot_df)
names(plot_df)[1] <- 'Lavi_True'
names(plot_df)[2] <- 'Lavi_Imputed'
plot_df <- na.omit(plot_df)
dim(plot_df)
plot_df_melt <- melt(as.matrix(plot_df))
head(plot_df_melt)
sum(plot_df_melt[,Var2=='Lavi_True'])
ggplot(plot_df_melt, aes(Var1, value)) +  geom_line(aes(colour = Var2)) + ylim(20,85) # endpoint_LAVI_MAX
ggplot(plot_df_melt, aes(Var1, value)) +  geom_line(aes(colour = Var2)) + ylim(16,2367) # endpoint_NTPROBNP
##### REGRESSIONS #####
