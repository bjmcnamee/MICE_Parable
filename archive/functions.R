clean_data <- function(data) {
        print('Cleaning data values...')
        # convert text to lower
        data <- data.frame(lapply(data, tolower))
        # remove leading or trailing white space
        data <- data.frame(apply(data, 2, function (x) gsub('^\\s+|\\s+$','',x)))
        # replace unknown labels with'NA'
        replacements <- c('Not done','#VALUE!','Unk','Not measured')
        for (r in replacements) {data <- data.frame(lapply(data, function(x) gsub(r, NA, x)))}
        # clean'N/A','Yes','No' values
        data <- data.frame(lapply(data, function(x) gsub('^n/a$|^na$|\\?', NA, x)))
        # AGAIN remove leading or trailing white space
        data <- data.frame(apply(data, 2, function (x) gsub('^\\s+|\\s+$','',x)))
        # convert text to numeric - no free text fields, saved as factor
        data <- as.data.frame(lapply(data, function(x) if(any(grepl('[[:digit:]]+\\.*[[:digit:]]*',x))) as.numeric(as.character(x)) else x))
        # convert NA text to true NA 
        make.true.NA <- function(x) {is.na(x) <- x==''; x}
        data <- data.frame(lapply(data, make.true.NA))
        print(paste('Dataset dimensions (rows x cols):', nrow(data),'x',ncol(data)))
        return (data) 
}

write_results <- function(name_var_impute, model, conf, Log_Transform, Exp_Transform){
        file <- paste(working_folder,'output/results.csv',sep='')
        write.table(name_var_impute, file, row.names=F, col.names=F, append=T, sep=',')
        write.table(data.frame('Regression Model', model), file, row.names=F, col.names=F, append=T, sep=',')
        write.table('Statistical Summary', file, row.names=F, col.names=F, append=T, sep=',')
        write.table(conf, file, row.names=F, col.names=T, append=T, sep=',')
        write.table('', file, row.names=F, col.names=F, append=T, sep=',')
        write.table('Log Transformation', file, row.names=F, col.names=F, append=T, sep=',')
        write.table(Log_Transform, file, row.names=F, col.names=T, append=T, sep=',')
        write.table('', file, row.names=F, col.names=F, append=T, sep=',')
        write.table('Exponential Transformation', file, row.names=F, col.names=F, append=T, sep=',')
        write.table(Exp_Transform, file, row.names=F, col.names=T, append=T, sep=',')
        write.table('', file, row.names=F, col.names=F, append=T, sep=',')
        write.table('', file, row.names=F, col.names=F, append=T, sep=',')
}
