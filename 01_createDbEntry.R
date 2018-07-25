# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 25/07/2018


## set variables and source libraries
source('header.R')

## connect to mysql database 
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and file table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe MetaFile;'))
cFileCol = dbGetQuery(db, paste('describe MetaFile;'))$Field[-1]

## load the metadata file
dfMeta = read.csv('dataExternal/LPMC metadata.csv', header=T, stringsAsFactors = F)
str(dfMeta)

## create the entry for samples
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=dfMeta$Sample.code, 
                       description= paste('group 1 is cytokine exposure', 'group2 is patient id', sep=';'),
                       group1 = dfMeta$Cytokine.exposure, group2= dfMeta$Patient.ID)
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# write this table to database
# dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 29;'))

cFileCol

dfCounts = read.csv('dataExternal/LPMC readcount.csv', header = T, row.names=1, sep='\t')
dim(dfCounts)

identical(colnames(dfCounts), as.character(dfSamples$title))
# put the names in the same order
table(colnames(dfCounts) %in% dfSamples$title)

i = match(dfSamples$title, colnames(dfCounts))
identical(dfSamples$title, colnames(dfCounts)[i])
dfCounts = dfCounts[,i]

# sanity check
identical(dfSamples$title, colnames(dfCounts))

n = make.names(paste('dfCounts for RNA-Seq data Jonathan data id 29 rds'))
n2 = paste0('~/Data/MetaData/', n)
# save(dfCounts, file=n2)
# 
# dbListFields(db, 'MetaFile')
# df = data.frame(idData=g_did, name=n, type='rds', location='~/Data/MetaData/',
#                 comment='data frame of count data from RNA-Seq experiment for Jonathan data id 29 rds')
# dbWriteTable(db, name = 'MetaFile', value=df, append=T, row.names=F)
# dbDisconnect(db)

