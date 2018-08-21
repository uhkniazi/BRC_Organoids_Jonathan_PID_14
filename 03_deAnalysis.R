# File: 03_deAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Modelling and selecting DE genes
# Date: 07/08/2018

## load the data
source('header.R')

## load the data
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 29) AND (MetaFile.comment like "%count%")')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
load(n)

## load the metadata i.e. covariates
q = paste0('select Sample.* from Sample where Sample.idData = 29')
dfSample = dbGetQuery(db, q)
dim(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)

## make count matrix
mCounts = as.matrix(dfCounts)
identical(colnames(mCounts), dfSample$title)

mData = mCounts
dim(mData)

# drop the samples where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
mData = mData[!(i< 3),]
dim(mData)
tail(rownames(mData))

## remove drop outs
ivProb = apply(mData, 1, function(inData) {
  inData[is.na(inData) | !is.finite(inData)] = 0
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})

hist(ivProb)
table(ivProb < 0.8)

mData = mData[!(ivProb < 0.8), ]
dim(mData)

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')
dim(mData.norm)
str(mData.norm)

## perform DE analysis
## delete sample section after testing
mData.norm = round(mData.norm, 0)
hist(mData.norm)

# set.seed(123)
# i = sample(1:nrow(mData.norm), 10, replace = F)
# dfData = data.frame(t(mData.norm[i,]))

dfData = data.frame(t(mData.norm))
dfData = stack(dfData)
dim(dfData)
dfData$fBatch = factor(dfSample$group1)
dfData$fAdjust1 = factor(dfSample$group2)
dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData$Coef.adj1 = factor(dfData$fAdjust1:dfData$ind)
## adjust the time as a continuous variable in each ind/gene
dim(dfData)
dfData = droplevels.data.frame(dfData)
dfData = dfData[order(dfData$Coef, dfData$Coef.adj1), ]
str(dfData)

# # # # # ## setup the model
# library(lme4)
# fit.lme1 = glmer.nb(values ~ 1  + (1 | Coef) + (1 | Coef.adj1), data=dfData)
# summary(fit.lme1)
# fit.lme2 = glmer.nb(values ~ 1 + (1 | Coef), data=dfData)
# summary(fit.lme2)
# 
# anova(fit.lme1, fit.lme2)
# 
# ran = ranef(fit.lme1, condVar=F)
# 
# plot(log(fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
# lines(lowess(log(fitted(fit.lme1)), resid(fit.lme1)), col=2)
# plot(density(log(dfData$values)))
# lines(density(log(fitted(fit.lme1))), col=2)
# lines(density(log(fitted(fit.lme2))), col=3)
# 
# plot(resid(fit.lme1), dfData$fBatch)
# plot(abs(resid(fit.lme1)), dfData$fBatch)
# plot(abs(resid(fit.lme1)), dfData$fAdjust1)
# plot(resid(fit.lme1), dfData$fAdjust1)

## setup the stan model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='nbinomResp2RandomEffects.stan')

## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(log(dfData$values+0.5)), 2*sd(log(dfData$values+0.5)))

# ### setup the input data for stan
# ## we give each gene its own variance term
# ## this however requires mapping of each variance term for each gene to the 
# ## number of coefficients/jitters for that gene
# l2 = split(dfData$Coef, dfData$ind)
# l2 = sapply(l2, function(x) length(unique(x)))
# l2 = gl(length(l2), unique(l2))

# lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
#                  Nclusters2=nlevels(dfData$Coef.adj1),
#                  Nclusters1_variance=nlevels(l2),
#                  Nsizes = nlevels(dfData$ind), ## each group of gene given its own dispersion term
#                  NgroupMap1=as.numeric(dfData$Coef),
#                  NgroupMap2=as.numeric(dfData$Coef.adj1),
#                  rGroupsJitter1Map=as.numeric(l2), ## mapping each gene variance term to its group of coefficients/jitters
#                  NsizeMap=as.numeric(dfData$ind),
#                  y=dfData$values, 
#                  gammaShape=l$shape, gammaRate=l$rate,
#                  intercept = mean(log(dfData$values+0.5)), intercept_sd= sd(log(dfData$values+0.5))*2)

lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 Nclusters2=nlevels(dfData$Coef.adj1),
                 NgroupMap1=as.numeric(dfData$Coef),
                 NgroupMap2=as.numeric(dfData$Coef.adj1),
                 Ncol=1, 
                 y=dfData$values, 
                 gammaShape=l$shape, gammaRate=l$rate,
                 intercept = mean(log(dfData$values+0.5)), intercept_sd= sd(log(dfData$values+0.5))*3)

fit.stan = sampling(stanDso, data=lStanData, iter=1500, chains=6,
                    pars=c('sigmaRan1', 'sigmaRan2', 'iSize',
                           'rGroupsJitter1'
                           ),
                    cores=6)
save(fit.stan, file='temp/fit.stan.nb_07Aug.rds')

print(fit.stan, c('sigmaRan1', 'sigmaRan2'), digits=3)
traceplot(fit.stan, 'sigmaRan1')
traceplot(fit.stan, 'sigmaRan2')
print(fit.stan, 'rGroupsJitter1')
print(fit.stan, 'iSize')

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)
# ## get the intercept at population level
#iIntercept = as.numeric(extract(fit.stan)$intercept)
## add the intercept to each random effect variable, to get the full coefficient
#mCoef = sweep(mCoef, 1, iIntercept, '+')

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

## split the data into the comparisons required
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$Coef))
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
head(d)
colnames(d) = c(colnames(d)[1:2], c('fBatch', 'ind'))
d$split = factor(d$ind)

levels(d$fBatch)
## repeat this for each comparison

## get a p-value for each comparison
l = tapply(d$cols, d$split, FUN = function(x, base='CD3', deflection='IL23+CD3') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif = getDifference(ivData = mCoef[,c[deflection]], ivBaseline = mCoef[,c[base]])
  r = data.frame(ind= as.character(d$ind[c[base]]), coef.base=mean(mCoef[,c[base]]), 
                 coef.deflection=mean(mCoef[,c[deflection]]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
  #return(format(r, digi=3))
  return(r)
})

dfResults = do.call(rbind, l)
dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

### plot the results
dfResults$logFC = dfResults$difference
dfResults$P.Value = dfResults$pvalue

library(org.Hs.eg.db)
df = select(org.Hs.eg.db, keys = as.character(rownames(dfResults)), columns = 'SYMBOL', keytype = 'ENSEMBL')
i = match(rownames(dfResults), df$ENSEMBL)
df = df[i,]
dfResults$SYMBOL = df$SYMBOL
identical(rownames(dfResults), df$ENSEMBL)
## produce the plots 
range(dfResults$difference)
f_plotVolcano(dfResults, 'IL23+CD3 vs CD3', fc.lim=c(-2.2, 2.2))

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfResults), names(m))
m = m[i]
identical(names(m), rownames(dfResults))
plotMeanFC(log(m), dfResults, 0.1, '')
table(dfResults$pvalue < 0.05)
table(dfResults$adj.P.Val < 0.01)
## save the results 
write.csv(dfResults, file='results/DEAnalysisIL23+CD3VsCD3.xls')



######### do a comparison with deseq2
f = factor(dfSample$group1)
f = relevel(f, 'Unstim')
dfDesign = data.frame(Treatment = f, fAdjust1 = factor(dfSample$group2), row.names=colnames(mData))

str(dfDesign)
oDseq = DESeqDataSetFromMatrix(mData, dfDesign, design = ~ Treatment + fAdjust1)
oDseq = DESeq(oDseq)
plotDispEsts(oDseq)
resultsNames(oDseq)
oRes = results(oDseq, contrast=c('Treatment', 'IL23.CD3', 'Unstim'))
plotMA(oRes)
temp = as.data.frame(oRes)
i = match(as.character(dfResults$ind), rownames(temp))
temp = temp[i,]
identical(as.character(dfResults$ind), rownames(temp))
plot(dfResults$logFC, temp$log2FoldChange, pch=20, main='ILC3: DEseq2 vs Bayes', xlab='Bayes', ylab='DEseq2')
table(oRes$padj < 0.01)
#write.csv(temp, file='results/DEAnalysisILC3VsSIO_Deseq2.xls')


#### plot the deseq2 volcano plot
dfResults = temp
dfResults$logFC = log(2^dfResults$log2FoldChange)
dfResults$P.Value = dfResults$pvalue
dfResults$adj.P.Val = dfResults$padj
head(rownames(dfResults))

df = select(org.Hs.eg.db, keys = as.character(rownames(dfResults)), columns = 'SYMBOL', keytype = 'ENSEMBL')
i = match(rownames(dfResults), df$ENSEMBL)
df = df[i,]
dfResults$SYMBOL = df$SYMBOL
identical(rownames(dfResults), df$ENSEMBL)## produce the plots 
f_plotVolcano(dfResults, 'DEseq2: ILC3 vs SIO', fc.lim=range(dfResults$logFC))
