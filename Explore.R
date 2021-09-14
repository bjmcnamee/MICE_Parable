##### EXPLORE DATA #####
library('gridExtra')
library("mice")
library(lattice)
working_folder <- '~/rstudio/scripts/Parable/' # edit this folder to match your working folder
source(paste(working_folder,'functions.R',sep=''))

# remove excess NAs, ie 10%+
data_lean <- drop_excess_NA_features(data, col_threshold=100, row_threshold=2.5, keep_endpoints=TRUE) # remove unwanted columns/observations
countNA <- sapply(data_lean, function(x) sum(is.na(x)))
c <- countNA[countNA>0]
names(data_lean)
length(c)
sum(c)
md.pat <- md.pattern(before_imputation_data, plot=TRUE, rotate.names=TRUE)


##### VISUAL INSPECTION #####
#imp1 <- imp_qp1; imp2 <- imp_qp2; imp3 <- imp_qp3

endpoints_LAVI <- unlist(features[1])
cd <- mice::complete(imp1)[, endpoints_LAVI]
mis <- is.na(imp1$data[,endpoints_LAVI[1]]) | is.na(imp1$data[,endpoints_LAVI[2]])
cd <- data.frame(mis = mis, cd)
x1 <- xyplot(jitter(cd[,endpoints_LAVI[1]], 10) ~ jitter(cd[,endpoints_LAVI[2]], 10),
             data = cd, groups = mis,
             main='1. full dataset',
             col = c(mdc(1), mdc(2)),
             xlab = "JDMRI_LAVI_MAX_BL",
             type = c("g","p"), ylab = "JDMRI_LAVI_MAX_18M",
             pch = c(1, 19),
             scales = list(alternating = 1, tck = c(1, 0)))
# '4. reduced dataset - less mice declared constant/collinear'
cd <- mice::complete(imp2)[, endpoints_LAVI]
mis <- is.na(imp2$data[,endpoints_LAVI[1]]) | is.na(imp2$data[,endpoints_LAVI[2]])
cd <- data.frame(mis = mis, cd)
x2 <- xyplot(jitter(cd[,endpoints_LAVI[1]], 10) ~ jitter(cd[,endpoints_LAVI[2]], 10),
             data = cd, groups = mis,
             main='2. reduced dataset',
             col = c(mdc(1), mdc(2)),
             xlab = "JDMRI_LAVI_MAX_BL",
             type = c("g","p"), ylab = "JDMRI_LAVI_MAX_18M",
             pch = c(1, 19),
             scales = list(alternating = 1, tck = c(1, 0)))
# '2. reduced dataset - less rows missing 10% values'
cd <- mice::complete(imp3)[, endpoints_LAVI]
mis <- is.na(imp3$data[,endpoints_LAVI[1]]) | is.na(imp3$data[,endpoints_LAVI[2]])
cd <- data.frame(mis = mis, cd)
x3 <- xyplot(jitter(cd[,endpoints_LAVI[1]], 10) ~ jitter(cd[,endpoints_LAVI[2]], 10),
             data = cd, groups = mis,
             main='3. reduced dataset',
             col = c(mdc(1), mdc(2)),
             xlab = "JDMRI_LAVI_MAX_BL",
             type = c("g","p"), ylab = "JDMRI_LAVI_MAX_18M",
             pch = c(1, 19),
             scales = list(alternating = 1, tck = c(1, 0)))
# '3. reduced dataset - less most extreme outliers'
cd <- mice::complete(imp4)[, endpoints_LAVI]
mis <- is.na(imp4$data[,endpoints_LAVI[1]]) | is.na(imp4$data[,endpoints_LAVI[2]])
cd <- data.frame(mis = mis, cd)
x4 <- xyplot(jitter(cd[,endpoints_LAVI[1]], 10) ~ jitter(cd[,endpoints_LAVI[2]], 10),
             data = cd, groups = mis,
             main='4. reduced dataset',
             col = c(mdc(1), mdc(2)),
             xlab = "JDMRI_LAVI_MAX_BL",
             type = c("g","p"), ylab = "JDMRI_LAVI_MAX_18M",
             pch = c(1, 19),
             scales = list(alternating = 1, tck = c(1, 0)))
grid.arrange(x1,x2,x3,x4,nrow = 2)

# endpoints_NTPROBNP
endpoints_NTPROBNP <- unlist(features[2])
imp <- imp2
NTPROBNP <- mice::complete(imp)[, endpoints_NTPROBNP]
for (i in 1:length(endpoints_NTPROBNP)) {
  mis <- is.na(imp$data[,endpoints_NTPROBNP[i]])
}
x1 <- xyplot(jitter(NTPROBNP[,2], 10) ~ jitter(NTPROBNP[,1], 10),
             data = cd, groups = mis,
             col = c(mdc(1), mdc(2)),
             xlab = substring(endpoints_NTPROBNP[2],2,),
             type = c("g","p"), ylab = endpoints_NTPROBNP[1],
             pch = c(1, 19), cex=0.5,
             scales = list(alternating = 1, tck = c(0, 0)))
x2 <- xyplot(jitter(NTPROBNP[,3], 10) ~ jitter(NTPROBNP[,1], 10),
             data = cd, groups = mis,
             col = c(mdc(1), mdc(2)),
             xlab = substring(endpoints_NTPROBNP[3],2,),
             type = c("g","p"), ylab = endpoints_NTPROBNP[1],
             pch = c(1, 19), cex=0.5,
             scales = list(alternating = 1, tck = c(0, 0)))
x3 <- xyplot(jitter(NTPROBNP[,4], 10) ~ jitter(NTPROBNP[,1], 10),
             data = cd, groups = mis,
             col = c(mdc(1), mdc(2)),
             xlab = substring(endpoints_NTPROBNP[4],2,),
             type = c("g","p"), ylab = endpoints_NTPROBNP[1],
             pch = c(1, 19), cex=0.5,
             scales = list(alternating = 1, tck = c(0, 0)))
x4 <- xyplot(jitter(NTPROBNP[,5], 10) ~ jitter(NTPROBNP[,1], 10),
             data = cd, groups = mis,
             col = c(mdc(1), mdc(2)),
             xlab = substring(endpoints_NTPROBNP[5],2,),
             type = c("g","p"), ylab = endpoints_NTPROBNP[1],
             pch = c(1, 19), cex=0.5,
             scales = list(alternating = 1, tck = c(0, 0)))
x5 <- xyplot(jitter(NTPROBNP[,6], 10) ~ jitter(NTPROBNP[,1], 10),
             data = cd, groups = mis,
             col = c(mdc(1), mdc(2)),
             xlab = substring(endpoints_NTPROBNP[6],2,),
             type = c("g","p"), ylab = endpoints_NTPROBNP[1],
             pch = c(1, 19), cex=0.5,
             scales = list(alternating = 1, tck = c(0, 0)))
x6 <- xyplot(jitter(NTPROBNP[,7], 10) ~ jitter(NTPROBNP[,1], 10),
             data = cd, groups = mis,
             col = c(mdc(1), mdc(2)),
             xlab = substring(endpoints_NTPROBNP[7],2,),
             type = c("g","p"), 
             ylab = endpoints_NTPROBNP[1],
             pch = c(1, 19), cex=0.5,
             scales = list(alternating = 1, tck = c(0, 0)))
grid.arrange(x1,x2,x3,x4,x5,x6,nrow = 2)

densityplot(imp1)

##### VISUAL INSPECTION #####
summary(data)

# list variables missing >0 NAs
countNA <- sapply(before_imputation_data, function(x) sum(is.na(x)))
countNA[countNA>0]

# list variables missing >0 NAs
countNA <- sapply(data, function(x) sum(is.na(x)))
c <- countNA[countNA>0]
names(c)
length(c)
sum(c)
length(data)
# remove excess NAs, ie 10%+
data_lean <- drop_excess_NA_features(data, col_threshold=10, row_threshold=10, keep_endpoints=TRUE) # remove unwanted columns/observations
countNA <- sapply(data_lean, function(x) sum(is.na(x)))
c <- countNA[countNA>0]
names(data_lean)
length(c)
sum(c)
md.pat <- md.pattern(before_imputation_data, plot=TRUE, rotate.names=TRUE)
# influx/outflux
fx <- flux(data, names(data))
plot(fx$influx, fx$outflux, xlim=c(0,1), ylim=c(0,1), abline(1,-1))
# how to add var name labels in position and legible https://stackoverflow.com/questions/7611169/intelligent-point-label-placement-in-r
row.names(fx)[fx$outflux<0.5] 
# pattern of missing data
md.pat <- md.pattern(data, plot=TRUE, rotate.names=TRUE)
write.csv(md.pat, '~/rstudio/scripts/Parable/output/md.pattern.before.csv', row.names=FALSE) # save cleaned dataframe
# structure
dim(data)
str(data, list.len=ncol(data))

# NA count
NA_count <- sum(is.na(data))
NA_not_count <- sum(!is.na(data))
NA_percent <- round(NA_count/(NA_count+NA_not_count)*100,1)
print(paste('NAs found ::', NA_count,'out of',NA_count+NA_not_count,'-', NA_percent,'% of dataset'))
# missing values
ini <- mice(data, maxit = 0)
missing_names <- names(ini$nmis[ini$nmis>0])
# barplot variable/column NAs by missing >10%
mp <- colSums(is.na(data))/(colSums(is.na(data))+colSums(!is.na(data)))*100 # variable missing %
mpt <- mp[mp > 5] # subset variable >10%
#mpt <- mpt[-40] # introduced in cleaning in error?
mpt <- mpt[order(names(mpt))] # get column names
mpt
length(mpt)
missing <- names(mpt)
missing
par(mar=c(5,11,4,1)+.1)
barplot(mpt, xlab = '% missing', main = 'Variables missing values > 5%', las=1, cex.names=0.7, horiz=TRUE)
# barplot observation/row NAs missing
NAs_by_row <- apply(data, 1, function(x) sum(is.na(x)))
barplot(NAs_by_row, xlim=c(0,list.len=ncol(data)), ylab = 'observations', xlab = 'missing count', main = 'Histogram of observations missing values', horiz=TRUE)

# compare list of variables where MICE applied pmm vs list with missing values
missing_names <- names(imp2$nmis[imp2$nmis>0])
missing_pmm <- names(imp2$method[imp2$method == 'pmm'])
setdiff(missing_names, missing_pmm)


library(VIM)
aggr(data, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
