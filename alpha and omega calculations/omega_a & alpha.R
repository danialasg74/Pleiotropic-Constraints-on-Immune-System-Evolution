library(stringr)
library(dunn.test)

# define Jukes-Cantor function
# from https://github.com/kousathanas/MultiDFE
# ***div.jukes***
# calculates Jukes-Cantor divergence.
# Input: x<-total sites,y<-site diffs
div.jukes<-function(x,y)
{
  d<-vector(length=length(x));
  for (i in 1:length(x))
  {
    if (y[i]<=0){d[i]=NA;next;}
    p=y[i]/x[i]
    if ((1-(4/3)*p)<0){d[i]=NA;next;}
    d[i]=(-3/4)*log(1-(4/3)*p)
  }
  
  return(d)
}



library(iMKT)
loadPopFly()
popral = subset(PopFlyData, Pop == "RAL")


# read in list of MultiDFE-outputted probabilities
setwd('.../MUN')


pleio_prob = read.table("pleio.txt")
colnames(pleio_prob) = c("replicate", "fixprob")
immun_prob = read.table("immun.txt")
colnames(immun_prob) = c("replicate", "fixprob")
devel_prob = read.table("devel.txt")
colnames(devel_prob) = c("replicate", "fixprob")
random_prob = read.table("random.txt")
colnames(random_prob) = c("replicate", "fixprob")



pleio_prob$class = "pleio"
immun_prob$class = "immun"
devel_prob$class = "devel"
random_prob$class = "random"
probs = rbind(pleio_prob, immun_prob, devel_prob, random_prob)

#directory containing bootstrap gene ids
setwd('.../bootstrap_genes')



probs$pi = 0 # from popfly
probs$p0 = 0 # from popfly
probs$di = 0 # from popfly
probs$d0 = 0 # from popfly
probs$mi = 0 # from popfly
probs$m0 = 0 # from popfly
probs$alpha = 0 # calculated
probs$omegaA = 0 # calculated
# add in columns for corrected data
probs$dicorr = 0
probs$d0corr = 0
probs$alphacorr = 0
probs$omegaAcorr = 0

rootdir = getwd()
prevclass = NULL

for (i in 1:nrow(probs)){
  bootstrapnum = probs[i,1] # bootstrap number
  geneidfile = paste(toString(bootstrapnum), ".csv", sep="")
  currentclass = sub("_.*", "", geneidfile) # go to the current class of gene (pleio, immun, devel)
  if (!identical(currentclass, prevclass)) {
    setwd(rootdir)            # go back to root
    setwd(currentclass)       # enter new class folder
    prevclass <- currentclass
  }
  geneinfo = read.table(geneidfile,header = TRUE) # read in file containing component genes  

  #initialize columns
  geneinfo$pi = 0 
  geneinfo$p0 = 0 
  geneinfo$di = 0 
  geneinfo$d0 = 0 
  geneinfo$mi = 0 
  geneinfo$m0 = 0 
  

  for (j in 1:nrow(geneinfo)){
    genename = geneinfo[j,1]
    rowmatch = match(genename, popral$Name)
    geneinfo[j,2] = popral[rowmatch, 6] 
    geneinfo[j,3] = popral[rowmatch, 5]
    geneinfo[j,4] = popral[rowmatch,7]
    geneinfo[j,5] = popral[rowmatch,8]
    geneinfo[j,6] = popral[rowmatch, 10]
    geneinfo[j,7] = popral[rowmatch, 11]
  }
  print(sum(!complete.cases(geneinfo)))
  #geneinfo = geneinfo[complete.cases(geneinfo), ]
  pi = sum(geneinfo$pi)
  p0 = sum(geneinfo$p0)
  di = sum(geneinfo$di)
  d0 = sum(geneinfo$d0)
  mi = sum(geneinfo$mi)
  m0 = sum(geneinfo$m0)  
  
  fixprob = probs[i, 2]
  probs[i,4] = pi
  probs[i,5] = p0
  probs[i,6] = di
  probs[i,7] = d0
  probs[i,8] = mi
  probs[i,9] = m0
  probs[i,10] = (di - (d0 * fixprob))/di # uncorrected alpha calculation (fixprob is w_na)
  probs[i,11] = (di - (d0 * fixprob))/d0 # uncorrected omega_a calculation
  dcorr = div.jukes(c(mi, m0), c(di, d0)) # calculate both dicorr and d0corr at once
  probs[i, 12] = dcorr[1] # add dicorr
  probs[i, 13] = dcorr[2] # add d0corr
  probs[i,14] = (dcorr[1] - (dcorr[2] * fixprob))/dcorr[1] # alpha calculation
  probs[i,15] = (dcorr[1] - (dcorr[2] * fixprob))/dcorr[2] # omega_a calculation
}




pleio_om_med = median(subset(probs, probs$class == 'pleio')$omegaAcorr)
immun_om_med = median(subset(probs, probs$class == 'immun')$omegaAcorr)
devel_om_med = median(subset(probs, probs$class == 'devel')$omegaAcorr)
rando_om_med = median(subset(probs, probs$class == 'random')$omegaAcorr)


###### plot omega_a

probs$class = factor(probs$class, levels = c('random','immun','pleio','devel'))
ggplot(probs, aes(x = class, y = omegaAcorr, fill = class)) +
  geom_violin(trim = TRUE, alpha = 1) +
  geom_boxplot(width = 0.1,outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  scale_fill_manual(
    values = c(
      immun = "#65c2a5",
      pleio = "#fc8d62",
      devel = "#8da0cb",
      random= 'grey'
    ),
    labels = c(
      immun = "Immunity",
      pleio = "Pleiotropic",
      devel = "Developmental",
      random = "Random"
      
    )
  )  +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "omega_a")+ylim(0,1)+scale_x_discrete(labels = c(
    immun = "Immunity",
    pleio = "Pleiotropic",
    devel = "Developmental",
    random = "Random"
  ))+annotate(
    "segment",
    x = 1.4, xend = 2.5,
    y = immun_om_med, yend = immun_om_med,
    colour = "black", # pick one of your fill colors
    linetype = "dashed",
    linewidth = 0.3
  )+annotate(
    "segment",
    x = 2.4, xend = 3.5,
    y = pleio_om_med, yend = pleio_om_med,
    colour = "black", # pick one of your fill colors
    linetype = "dashed",
    linewidth = 0.3
  )+annotate(
    "segment",
    x = 3.3, xend = 4.4,
    y = devel_om_med, yend = devel_om_med,
    colour = "black", # pick one of your fill colors
    linetype = "dashed",
    linewidth = 0.3
  )+annotate(
    "segment",
    x = 0.5, xend = 1.6,
    y = rando_om_med, yend = rando_om_med,
    colour = "black", # pick one of your fill colors
    linetype = "dashed",
    linewidth = 0.3
  )+ annotate(
    "text",
    x = 1,
    y  = 0.5,
    label = "c",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 2,
    y  = 0.5,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 3,
    y = 0.5,
    label = "ab",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 4,
    y = 0.5,
    label = "b",
    colour = "black",
    size = 5)


####### stat test
x = c(probs$omegaAcorr)
grp = c(probs$class)
kruskal.test(x ~ grp)
pairwise.wilcox.test(x, grp, p.adjust.method = "BH", exact = FALSE, na.action = na.omit)

#random = c  
#immun  = a  
#pleio  = ab  
#devel  = b  
 




















###### plot alpha

pleio_al_med = median(subset(probs, probs$class == 'pleio')$alphacorr)
immun_al_med = median(subset(probs, probs$class == 'immun')$alphacorr)
devel_al_med = median(subset(probs, probs$class == 'devel')$alphacorr)
rando_al_med = median(subset(probs, probs$class == 'random')$alphacorr)

ggplot(probs, aes(x = class, y = alphacorr, fill = class)) +
  geom_violin(trim = TRUE, alpha = 1) +
  geom_boxplot(width = 0.1,outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  scale_fill_manual(
    values = c(
      immun = "#65c2a5",
      pleio = "#fc8d62",
      devel = "#8da0cb",
      random = 'grey'
    ),
    labels = c(
      immun = "Immunity",
      pleio = "Pleiotropic",
      devel = "Developmental",
      random = "Random"
    )
  )  +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "alpha")+ylim(0,1)+scale_x_discrete(labels = c(
    immun = "Immunity",
    pleio = "Pleiotropic",
    devel = "Developmental",
    random = "Random"
  ))+annotate(
    "segment",
    x = 1.4, xend = 2.5,
    y = immun_al_med, yend = immun_al_med,
    colour = "black", # pick one of your fill colors
    linetype = "dashed",
    linewidth = 0.3
  )+annotate(
    "segment",
    x = 2.4, xend = 3.5,
    y = pleio_al_med, yend = pleio_al_med,
    colour = "black", # pick one of your fill colors
    linetype = "dashed",
    linewidth = 0.3
  )+annotate(
    "segment",
    x = 3.3, xend = 4.4,
    y = devel_al_med, yend = devel_al_med,
    colour = "black", # pick one of your fill colors
    linetype = "dashed",
    linewidth = 0.3
  )+annotate(
    "segment",
    x = 0.5, xend = 1.6,
    y = rando_al_med, yend = rando_al_med,
    colour = "black", # pick one of your fill colors
    linetype = "dashed",
    linewidth = 0.3
  )+ annotate(
    "text",
    x = 1,
    y  = 0.5,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 2,
    y  = 0.5,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 3,
    y = 0.5,
    label = "b",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 4,
    y = 0.5,
    label = "c",
    colour = "black",
    size = 5)

####### stat test
x = c(probs$alphacorr)
grp = c(probs$class)
kruskal.test(x ~ grp)
pairwise.wilcox.test(x, grp, p.adjust.method = "BH", exact = FALSE, na.action = na.omit)


#random = a  
#immun  = a  
#pleio  = b  
#devel  = c  
