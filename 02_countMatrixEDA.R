# File: 02_countMatrixEDA.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: quality checks on the count matrix with covariates
# Date: 25/07/2018

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

## load CDiagnostics and test
## compare the normalised and raw data
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(log(mData.norm+1), 'Normalised')
oDiag.2 = CDiagnosticPlots(log(mData+1), 'Original')

# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
str(dfSample)
fBatch = factor(dfSample$group1)
fBatch = factor(dfSample$group2)

#pdf('results/matrixClustering.pdf')

## compare the 2 methods using various plots
par(mfrow=c(1,2))
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)
plot.missing.summary(oDiag.2, fBatch, axis.label.cex = 0.7, cex.main=1)

plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.PCA(oDiag.2, fBatch, cex.main=1)
# 
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.8, cex.main=0.7)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.1 = CDiagnosticPlotsSetParameters(oDiag.1, l)
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)
plot.PCA(oDiag.1, fBatch)
plot.PCA(oDiag.2, fBatch)
plot.PCA(oDiag.1, fBatch, csLabels = dfSample$group2)
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.7)

# ## extreme values
# oDiag.1 = calculateExtremeValues(oDiag.1)
# oDiag.2 = calculateExtremeValues(oDiag.2)
# m1 = mGetExtremeValues(oDiag.1)
# m2 = mGetExtremeValues(oDiag.2)
# 
# ## samples with most extreme values
# apply(m1, 2, function(x) sum(x > 0))
# apply(m2, 2, function(x) sum(x > 0))
# 
# ## variables that are contributing to this
# v1 = apply(m1, 1, function(x) sum(x > 0))
# v2 = apply(m2, 1, function(x) sum(x > 0))
# 
# which(v1 > 0)
# which(v2 > 0)

###############################################################################################
################ clusters and covariates explaining variance
mPC = oDiag.1@lData$PCA$x[,1:2]
i = which(mPC[,2] > 50)
fNewBatch = rep(1, times=length(dfSample$title))
fNewBatch[i] = 2
fNewBatch = factor(fNewBatch)

## check if batch assigned correctly
plot.PCA(oDiag.1, fNewBatch)

## try a linear mixed effect model to account for varince
library(lme4)

dfData = data.frame(mPC)
dfData = stack(dfData)
dfData$fBatch = factor(dfSample$group1)
dfData$fAdjust1 = factor(dfSample$group2)
dfData$fAdjust2 = fNewBatch

dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData$Coef.adj1 = factor(dfData$fAdjust1:dfData$ind)
dfData$Coef.adj2 = factor(dfData$fAdjust2:dfData$ind)

str(dfData)

fit.lme1 = lmer(values ~ 1  + (1 | Coef) + (1 | Coef.adj1) + (1 | Coef.adj2), data=dfData)
summary(fit.lme1)

fit.lme2 = lmer(values ~ 1  + (1 | Coef) + (1 | Coef.adj1), data=dfData)
summary(fit.lme2)

fit.lme3 = lmer(values ~ 1  + (1 | Coef), data=dfData)
summary(fit.lme3)

anova(fit.lme1, fit.lme2, fit.lme3)

par(p.old)
plot((fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
lines(lowess((fitted(fit.lme1)), resid(fit.lme1)), col=2)
hist(dfData$values, prob=T)
lines(density(fitted(fit.lme1)))

plot((fitted(fit.lme2)), resid(fit.lme2), pch=20, cex=0.7)
lines(lowess((fitted(fit.lme2)), resid(fit.lme2)), col=2)
hist(dfData$values, prob=T)
lines(density(fitted(fit.lme2)))
lines(density(fitted(fit.lme1)), col=2)

plot((fitted(fit.lme3)), resid(fit.lme3), pch=20, cex=0.7)
lines(lowess((fitted(fit.lme3)), resid(fit.lme3)), col=2)
hist(dfData$values, prob=T)
lines(density(fitted(fit.lme2)))
lines(density(fitted(fit.lme1)), col=2)
lines(density(fitted(fit.lme3)), col=3)

## fit model with stan
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='tResponse2RandomEffectsNoFixed.stan')


### try a t model without mixture
lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 Nclusters2=nlevels(dfData$Coef.adj1),
                 NgroupMap1=as.numeric(dfData$Coef),
                 NgroupMap2=as.numeric(dfData$Coef.adj1),
                 y=dfData$values, 
                 gammaShape=0.5, gammaRate=1e-4)

fit.stan = sampling(stanDso, data=lStanData, iter=10000, chains=4,
                    pars=c('betas', 'sigmaRan1', 'sigmaRan2',
                           'nu', 'sigmaPop', 'mu',
                           'rGroupsJitter1', 'rGroupsJitter2'),
                    cores=4)#, init=initf, control=list(adapt_delta=0.99, max_treedepth = 12))

print(fit.stan, c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaPop', 'nu'), digits=3)
traceplot(fit.stan, 'betas')
traceplot(fit.stan, 'sigmaRan1')
traceplot(fit.stan, 'sigmaRan2')


m = cbind(extract(fit.stan)$sigmaRan1, extract(fit.stan)$sigmaRan2) 
dim(m)
colnames(m) = c('Treatment', 'Patient')
pairs(m, pch=20, col='lightgrey', cex=0.5)
m = log(m)

library(lattice)
df = stack(data.frame(m))
histogram(~ values | ind, data=df, xlab='Log SD')

hist(dfData$values, prob=T)
mFitted = extract(fit.stan)$mu

apply(mFitted[sample(1:nrow(mFitted), size = 100), ], 1, function(x) lines(density(x)))

