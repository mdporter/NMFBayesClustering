###############################################################################
# Analysis of crabs data (Section 5.3) in 
#   "Optimal Bayesian Clustering using Non-negative Matrix Factorization"
#   by Ketong Wang and Michael D. Porter
# This code produces Figure 8 and Table 4
# April 30, 2018
###############################################################################

#-- Required Packages
library(mcclust)
library(mcclust.ext)
library(NMF)
library(MASS)
library(dplyr)


#==============================================================#
#-- Get Data
#==============================================================#
# Note: must have crabs-draws.rds in path

#-- crab data  
data(crabs, package="MASS")   # load crabs data
true.cl = paste(crabs$sp, crabs$sex, sep=':')   # true classes

#-- get cluster (hard-partition) data
crabs.draw = readRDS('crabs-draws.rds')

#-- make pairwise similarity matrix (\pi)
psm = mcclust::comp.psm(crabs.draw)
n = nrow(psm)

#==============================================================#
#-- Get optimal partitions from existing methods
#==============================================================#
library(mcclust)
library(mcclust.ext)

max.k = 10    # maximum number of clusters

#- MinBinder
mbind = mcclust::minbinder(psm, max.k = max.k, method="avg")

#- MaxPear
mpear = mcclust::maxpear(psm, max.k = max.k, method="avg")

#- minVI
min_VI = mcclust.ext::minVI(psm, max.k = max.k, method="avg")

#- Medvedovic method
medvedovic = matrix( mcclust::medv(psm, h=0.99), nrow=1 )
medvedovic$cl = as.vector(medvedovic)


#==============================================================#
#-- NMF
#==============================================================#
library(NMF)
library(tibble)
  
#-- NMF settings
# note: nmf() uses parallel unless set .pbackend = NA
ncores = NA

K_seq = 2:max.k
NMF.seed = 123456
NMF.iters = 100

#-- Run all NMF models
algs_list <- c('brunet', 'lee', 'nsNMF', 'offset')
NMF_results = tibble()

for(algName in algs_list) {
  for(K in K_seq) {
    cat('NMF algorithm:', algName,  '| K = ', K, '\n')    
    myAlg <- nmfAlgorithm(algName)
    res <- nmf(psm, rank=K, method=myAlg, seed=NMF.seed, nrun=NMF.iters, .pbackend=ncores)
    H = coef(res)
    cls = apply(H, 2, which.max)
    NMF_results = bind_rows(NMF_results,
                            tibble(obs=1:n, cls, K, method=algName))
  }
}    

NMF_results = NMF_results %>% 
  mutate(method = recode(method, lee = 'NMF-ls', brunet='NMF-kl', 
                                                     nsNMF='NMF-ns', offset='NMF-offset'))



#- Evaluate over rank K
NMF = NMF_results %>%
  group_by(method, K) %>%
  do(data.frame(score=mcclust::pear(.$cls, psm), 
                K.actual = n_distinct(.$cls))) %>% 
  group_by(method) %>% filter(score == max(score)) %>% 
  filter(K==min(K))

#- Get optimal partitions
NMF_cl = NMF_results %>% 
  semi_join(NMF, by=c("method", "K")) %>% 
  split(f=.$method)

#- re-label clusters
NMF_cl = NMF_cl %>% 
  lapply(function(x)  x %>% arrange(obs) %>% mutate(cl = factor(cls, levels=unique(cls)),
                                                    cl = as.integer(cl)))



#==============================================================#
#-- Get partitions from all methods
#==============================================================#

CLS = bind_cols(
  sapply(NMF_cl, '[[', 'cl') %>% as.data.frame , 
  data.frame(MinBinder = mbind$cl, 
             MaxPEAR = mpear$cl, 
             MinVI = min_VI$cl, 
             Medv = medvedovic$cl)
) %>% as_tibble() %>%  
  rename_all(function(x) recode(x, 
                                lee = 'NMF-ls', brunet='NMF-kl', 
                                nsNMF='NMF-ns', offset='NMF-offset')) %>% 
  select(`NMF-ls`, everything()) %>%    # change order of columns
  as_tibble()



#==============================================================#
#-- Figure 8 from paper: Partition Matrices (all methods)
#==============================================================#

#-- Plot for Crab data
plot_crabs <- function(psm, ...){
  image(1:n, 1:n, 1-psm, col = heat.colors(20),
        xaxt='n', yaxt='n', xlab='', ylab='', 
        cex.main=2, ...)
  axis(1, at=seq(0, n, by=50), cex.axis=1.5)
  axis(2, at=seq(0, n, by=50), cex.axis=1.5, las=1)
  box(col="lightgray")  
  for(d in c(0,50,100,150)){
    polygon(x=d+c(0,0, 50, 50)+.5, y=d+c(0,50, 50, 0)+.5, lty=2, lwd=3, xpd=TRUE)
  }
}


#-- Function to convert partition to similarity matrix
cl2mat <- function(cl) outer(cl, cl, function(x,y) 1*(x==y))



#-- Find if any partitions are identical

data.frame(
  pair = combn(colnames(CLS), 2) %>% apply(2, paste, collapse=":"),
  ndiff = combn(CLS, 2, FUN=function(x) sum(x[,1]!=x[,2]))
)

# NMF-ls:NMF-ns:Medv 
# NMF-kl:MaxPEAR:MinVI


#-- Figure 8  

layout(matrix(1:8, ncol=4, nrow=2), respect=TRUE)
par(mar=c(4, 4, 4, 2))

plot_crabs(psm, main='MCMC')
plot_crabs(cl2mat(CLS$MinBinder), 
           main=paste0("MinBinder (K=", n_distinct(CLS$MinBinder), ")"))
plot_crabs(cl2mat(CLS$MaxPEAR), 
           main=paste0("MaxPear/MinVI (K=", n_distinct(CLS$MaxPEAR), ")"))
# plot_crabs(cl2mat(CLS$MinVI), main=paste0("MinVI (K=", n_distinct(CLS$MinVI), ")"))
plot_crabs(cl2mat(CLS$Medv), 
           main=paste0("Med (K=", n_distinct(CLS$Medv), ")"))
plot_crabs(cl2mat(CLS$'NMF-ls'),
           main=paste0("NMF-ls (K=", n_distinct(CLS$'NMF-ls'), ")"))
plot_crabs(cl2mat(CLS$'NMF-ns'),
           main=paste0("NMF-ns (K=", n_distinct(CLS$'NMF-ns'), ")"))
plot_crabs(cl2mat(CLS$'NMF-kl'),
           main=paste0("NMF-kl (K=", n_distinct(CLS$'NMF-kl'), ")"))
plot_crabs(cl2mat(CLS$'NMF-offset'),
           main=paste0("NMF-offset (K=", n_distinct(CLS$'NMF-offset'), ")"))


#==============================================================#
#-- Table 4 from paper: performance table (Rand, AR, VI)
#==============================================================#

#-- Make performance table
apply(CLS, 2, function(cls) 
  c(Rand = mcclust::arandi(cls, cl2=true.cl, adjust=FALSE),
    AR = mcclust::arandi(cls, cl2=true.cl, adjust=TRUE),
    VI = mcclust::vi.dist(cls, cl2=true.cl))) %>% 
  round(3)


  









