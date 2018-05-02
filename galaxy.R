###############################################################################
# Analysis of galaxies data (Section 5.2) in 
#   "Optimal Bayesian Clustering using Non-negative Matrix Factorization"
#   by Ketong Wang and Michael D. Porter
# This code produces Figures 6 and 7 
# April 30, 2018
###############################################################################

#-- Required Packages

library(tidyverse)
library(mcclust.ext)
library(mcclust)
library(NMF)

#==============================================================#
#-- Get Data
#==============================================================#

#-- Load Galaxy Data
data(galaxy.draw, package="mcclust.ext")
data(galaxy.fit, package="mcclust.ext")


#-- Get pairwise similarity matrix (\pi)
psm = mcclust::comp.psm(galaxy.draw)
n = nrow(psm)

image(1:n, 1:n, 1-psm, col = heat.colors(20),
      xaxt='n', yaxt='n', xlab='', ylab='')
axis(1, at=seq(0, n, by=5))
axis(2, at=seq(0, n, by=5), las=1)
abline(v=0:n-.5, col="grey90", lwd=.5)
abline(h=0:n-.5, col="grey90", lwd=.5)
box()




#==============================================================#
#-- Get optimal partitions from existing methods
#==============================================================#

max.k = 16    # maximum number of clusters

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
#-- Get optimal partitions from NMF methods
#==============================================================#

#- NMF settings
# note: nmf() uses parallel unless set .pbackend = NA
ncores = NA

K_seq = 2:max.k
NMF.seed = 123456
NMF.iters = 100

#- Run all NMF models
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
  do(data.frame(score=mcclust::pear(.$cls, psm))) %>% 
  group_by(method) %>% filter(score == max(score)) 


#-- Get optimal partitions
NMF_cl = NMF_results %>% 
  semi_join(NMF, by=c("method", "K")) %>% 
  split(f=.$method)


#- re-label clusters
NMF_cl = NMF_cl %>% 
  lapply(function(x)  x %>% arrange(obs) %>% mutate(cl = factor(cls, levels=unique(cls)),
                                                    cl = as.integer(cl)))




#==============================================================#
#-- Figure 6 from Paper: Clustering Solutions
#==============================================================#

#-- Settings
x = galaxy.fit$x   
y = rep(0, length(x))
symbols <-  c(16, 0, 17, 3, 15, 4, 1)
colors <- c("dodgerblue2", "red3", "green3", "slateblue", "darkorange", 
            "skyblue1", "violetred4") 



#-- Plot
par(mfrow=c(2,4), mar=c(3,1,3,1))

with(mbind, plot(x,y, pch=symbols[cl], col=colors[cl], 
                 xlab='', ylab='', yaxt='n', 
                 main=paste0('MinBinder (K=', n_distinct(cl),')')))

with(mpear, plot(x,y, pch=symbols[cl], col=colors[cl],
                 xlab='', ylab='', yaxt='n',
                 main=paste0('MaxPEAR (K=', n_distinct(cl),')')))

with(medvedovic, plot(x,y, pch=symbols[cl], col=colors[cl],
                      xlab='', ylab='', yaxt='n',
                      main=paste0('Medvedovic (K=', n_distinct(cl),')')))

with(NMF_cl$'NMF-ns', 
     plot(x,y, pch=symbols[cl], col=colors[cl], 
          xlab='', ylab='', yaxt='n',                   
          main=paste0('NMF-ns (K=', n_distinct(cl),')')))

with(min_VI, plot(x,y, pch=symbols[cl], col=colors[cl], 
                  xlab='', ylab='', yaxt='n',                   
                  main=paste0('MinVI (K=', n_distinct(cl),')')))

with(NMF_cl$'NMF-kl', 
     plot(x,y, pch=symbols[cl], col=colors[cl], 
          xlab='', ylab='', yaxt='n',                   
          main=paste0('NMF-kl (K=', n_distinct(cl),')')))

with(NMF_cl$'NMF-ls', 
     plot(x,y, pch=symbols[cl], col=colors[cl], 
          xlab='', ylab='', yaxt='n',                   
          main=paste0('NMF-ls (K=', n_distinct(cl),')')))

with(NMF_cl$'NMF-offset', 
     plot(x,y, pch=symbols[cl], col=colors[cl], 
          xlab='', ylab='', yaxt='n',                   
          main=paste0('NFM-offset (K=', n_distinct(cl),')')))


#==============================================================#
#-- Figure 7 from Paper: Soft-clustering with NMF-kl 
#==============================================================#

#- NMF settings
# use same values as above

#-- Run brunet (NMF-kl) 
algName = 'brunet'
K_seq = 3:5

S = vector('list', length(K_seq)); names(S) = paste0("K=", K_seq)
for(i in seq_along(K_seq)) {
  K = K_seq[i]
  print(paste('starting K =', K))
  myAlg <- nmfAlgorithm(algName)
  res <- nmf(psm, rank=K, method=myAlg, seed=NMF.seed, nrun=NMF.iters, .pbackend=ncores)
  H = coef(res); colnames(H) = 1:ncol(H)
  r = rle(apply(H, 2, which.max))    
  ord = r$values[1:nrow(H)]
  H = H[ord,]          # order clusters by x values
  S[[i]] = apply(H, 2, function(x) x/sum(x))
}


#-- Set observations to label
xlabs = rep("", n)
xlabs[c(5, 8, 15, 25, 35, 45, 55, 65, 75, 78, 82 )] =c(5, 8, 15, 25, 35, 45, 55, 65, 75, 78, 82 )

#-- Plot
bind_rows(  
  t(S[["K=3"]]) %>% as_tibble() %>% 
    mutate(obs=row_number()) %>% 
    gather(cluster, prob, -obs) %>% 
    mutate(cluster=stringr::str_sub(cluster, 2 )) %>% 
    add_column(K=paste0("K=",3)), 
  t(S[["K=4"]]) %>% as_tibble() %>% 
    mutate(obs=row_number()) %>% 
    gather(cluster, prob, -obs) %>% 
    mutate(cluster=stringr::str_sub(cluster, 2 )) %>% 
    add_column(K=paste0("K=",4)),
  t(S[["K=5"]]) %>% as_tibble() %>% 
    mutate(obs=row_number()) %>% 
    gather(cluster, prob, -obs) %>% 
    mutate(cluster=stringr::str_sub(cluster, 2 )) %>% 
    add_column(K=paste0("K=",5))
) %>% 
  ggplot(aes(factor(obs), prob)) + 
  geom_col(aes(fill=cluster, group=factor(cluster, sort(unique(cluster), decreasing=TRUE))), color="black")+ 
  geom_hline(yintercept=c(.20, .40, .6, .80), color='grey80')  + 
  scale_y_continuous(expand = c(0, 0), breaks=c(.20, .40, .6, .80)) + 
  scale_fill_manual(values=colors)   + 
  scale_x_discrete(labels=xlabs) +
  facet_wrap(~K, ncol=1) + 
  theme(
    plot.title = element_text(hjust = 0.5),  
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank())


