library(jsonlite)
library(dplyr)
library(stringr)
library(ggplot2)


################################### Processing dN/dS for the 4 classes of genes ###################################

###################################################################################################################


##### file for anntation ####

ncbi_ann = read.delim('.../ncbi_droso_ann.tsv')
ncbi_ann = ncbi_ann[,c(6,7,12)]
names(ncbi_ann)[3] = 'protein'








#################### pleiotropic genes ####################

setwd('.../pleio')

# read JSON file
pleio = fromJSON("dN_and_ds.txt")

# convert the nested list to a flat data frame
pleio = do.call(rbind, lapply(names(pleio), function(name) {
  cbind(file = name, as.data.frame(pleio[[name]]))
}))
rownames(pleio) = NULL

# extract protein and domain names using regex
library(stringr)

pleio$protein = str_extract(pleio$file, "NP_[0-9]+\\.[0-9]+")
pleio$domain  = str_extract(pleio$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")

# for rows without a domain (just the main protein), set domain = "none"
pleio$domain[is.na(pleio$domain)] <- "none"

# reorder columns
pleio = pleio[, c("protein", "domain", "dN", "dS")]


pleio = merge(pleio, ncbi_ann, by = 'protein')
pleio$W = (pleio$dN)/(pleio$dS)
pleio$kind  = ifelse(pleio$domain == "none", "Gene", "Domain")
pleio$kind = factor(pleio$kind, levels = c('Gene' , 'Domain'))








#################### immune genes ####################


setwd('.../immun')


# read JSON file
immun = fromJSON("dN_and_ds.txt")

# convert the nested list to a flat data frame
immun = do.call(rbind, lapply(names(immun), function(name) {
  cbind(file = name, as.data.frame(immun[[name]]))
}))
rownames(immun) = NULL

# extract protein and domain names using regex
library(stringr)

immun$protein = str_extract(immun$file, "NP_[0-9]+\\.[0-9]+")
immun$domain  = str_extract(immun$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")

# for rows without a domain (just the main protein), set domain = "none"
immun$domain[is.na(immun$domain)] <- "none"

# reorder columns
immun = immun[, c("protein", "domain", "dN", "dS")]




immun = merge(immun, ncbi_ann, by = 'protein')
immun$W = (immun$dN)/(immun$dS)
immun$kind  = ifelse(immun$domain == "none", "Gene", "Domain")
immun$kind = factor(immun$kind, levels = c('Gene' , 'Domain'))








#################### developmental genes ####################


setwd('.../devel')


# read JSON file
devel = fromJSON("dN_and_ds.txt")

# convert the nested list to a flat data frame
devel = do.call(rbind, lapply(names(devel), function(name) {
  cbind(file = name, as.data.frame(devel[[name]]))
}))
rownames(devel) = NULL

# extract protein and domain names using regex
library(stringr)

devel$protein = str_extract(devel$file, "NP_[0-9]+\\.[0-9]+")
devel$domain  = str_extract(devel$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")

# for rows without a domain (just the main protein), set domain = "none"
devel$domain[is.na(devel$domain)] <- "none"

# reorder columns
devel = devel[, c("protein", "domain", "dN", "dS")]




devel = merge(devel, ncbi_ann, by = 'protein')
devel$W = (devel$dN)/(devel$dS)
devel$kind  = ifelse(devel$domain == "none", "Gene", "Domain")
devel$kind = factor(devel$kind, levels = c('Gene' , 'Domain'))








#################### random genes ####################


setwd('.../Random')


# read JSON file
rando = fromJSON("dN_and_ds.txt")

# convert the nested list to a flat data frame
rando = do.call(rbind, lapply(names(rando), function(name) {
  cbind(file = name, as.data.frame(rando[[name]]))
}))
rownames(rando) = NULL

# extract protein and domain names using regex
library(stringr)

rando$protein = str_extract(rando$file, "NP_[0-9]+\\.[0-9]+")
rando$domain  = str_extract(rando$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")

# for rows without a domain (just the main protein), set domain = "none"
rando$domain[is.na(rando$domain)] <- "none"

# reorder columns
rando = rando[, c("protein", "domain", "dN", "dS")]




rando = merge(rando, ncbi_ann, by = 'protein')
rando$W = (rando$dN)/(rando$dS)
rando$kind  = ifelse(rando$domain == "none", "Gene", "Domain")
rando$kind = factor(rando$kind, levels = c('Gene' , 'Domain'))








##################### plot dN/dS (domain = y-axis; gene = x-axis) #################################################

###################################################################################################################


############ combine and filter out ds <3 ########

pleio$group = 'Pleiotropic'
immun$group = 'Immune'
devel$group = 'Developmental'
rando$group = 'Random'

All = rbind(pleio,immun,devel,rando)
All_uns = subset(All, All$dS < 3)


rm(pleio, immun, devel,rando)


############ Add gene-level dN/dS to domain data ########

All_uns$W_gene = NA

for (i in seq_len(nrow(All_uns))) {
  PROT = subset(All_uns, All_uns$protein == All_uns$protein[i])
  
  w_val <- subset(PROT, PROT$domain == "none")$W
  
  if (length(w_val) > 0) {
    # if there are multiple, just take the first; adjust if you want something else
    All_uns$W_gene[All_uns$protein == All_uns$protein[i]] <- w_val[1]
  } else {
    cat("Skipped index:", i, " protein:", All_uns$protein[i], "\n")
    next
  }
}


All_Dom_uns = subset(All_uns, All_uns$kind=='Domain')
names(All_Dom_uns)[7] = "W_domain"

############ Break the dataframe into 4 groups ########

Pleiotropic_uns = subset(All_Dom_uns,  All_Dom_uns$group == 'Pleiotropic')
Developmental_uns = subset(All_Dom_uns,  All_Dom_uns$group == 'Developmental')
Immune_uns = subset(All_Dom_uns,  All_Dom_uns$group == 'Immune')
Random_uns = subset(All_Dom_uns,  All_Dom_uns$group == 'Random')


Pleiotropic_uns = na.omit(Pleiotropic_uns)
Immune_uns = na.omit(Immune_uns)
Developmental_uns = na.omit(Developmental_uns)
Random_uns = na.omit(Random_uns)

############ Calculate % with Y[i] > X[i] ############ 

above_line_pleio = Pleiotropic_uns[Pleiotropic_uns$W_domain > Pleiotropic_uns$W_gene, ]
above_line_immun = Immune_uns[Immune_uns$W_domain > Immune_uns$W_gene, ]
above_line_devel = Developmental_uns[Developmental_uns$W_domain > Developmental_uns$W_gene, ]
above_line_rando = Random_uns[Random_uns$W_domain > Random_uns$W_gene, ]


pleio_frac = nrow(above_line_pleio)/nrow(Pleiotropic_uns)
immun_frac = nrow(above_line_immun)/nrow(Immune_uns)
devel_frac = nrow(above_line_devel)/nrow(Developmental_uns)
rando_frac = nrow(above_line_rando)/nrow(Random_uns)



######################## Plot Fig1 A ########################

pleio_plot_uns = ggplot(Pleiotropic_uns, aes(x = W_gene, y = W_domain)) +
  geom_point(size = 0.5, alpha = 1, color = "#fc8d62") +
  geom_smooth(method = "lm", color = "#fc8d62") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.2), limits = c(0, 0.6)) +
  xlab("dN/dS Gene") +
  ylab("") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    title = element_text(size = 14)
  )+ggtitle("Pleiotropic")+ annotate(
    "text",
    x = 0,
    y = 0.6,
    label = paste0("% with Y[i] > X[i]: ", round(100 * pleio_frac, 1), "%"),
    
    hjust = 0,
    vjust = 1,
    size = 5
  )

immun_plot_uns = ggplot(Immune_uns, aes(x = W_gene, y = W_domain)) +
  geom_point(size = 0.5, alpha = 1, color = "#65c2a5") +
  geom_smooth(method = "lm", color = "#65c2a5") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.2), limits = c(0, 0.6)) +
  xlab("dN/dS Gene") +
  ylab("") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    title = element_text(size = 14)
  )+ggtitle("Immunity")+ annotate(
    "text",
    x = 0,
    y = 0.6,
    label = paste0("% with Y[i] > X[i]: ", round(100 * immun_frac, 1), "%"),
    
    hjust = 0,
    vjust = 1,
    size = 5
  )

devel_plot_uns = ggplot(Developmental_uns, aes(x = W_gene, y = W_domain)) +
  geom_point(size = 0.5, alpha = 1, color = "#8da0cb") +
  geom_smooth(method = "lm", color = "#8da0cb") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.2), limits = c(0, 0.6)) +
  xlab("dN/dS Gene") +
  ylab("") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    title = element_text(size = 14)
  )+ggtitle("Developmental")+ylim(0,0.6)+ annotate(
    "text",
    x = 0,
    y = 0.6,
    label = paste0("% with Y[i] > X[i]: ", round(100 * devel_frac, 1), "%"),
    
    hjust = 0,
    vjust = 1,
    size = 5
  )



rando_plot_uns = ggplot(Random_uns, aes(x = W_gene, y = W_domain)) +
  geom_point(size = 0.5, alpha = 1, color = "grey") +
  geom_smooth(method = "lm", color = "grey") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  annotate("text", x = 0.3, y = 0.45, label = "y = x", color = "black", size = 8, fontface = "italic") +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.2), limits = c(0, 0.6)) +
  xlab("dN/dS Gene") +
  ylab("dN/dS Domain") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    title = element_text(size = 14)
  )+ggtitle("Random")+ annotate(
    "text",
    x = 0,
    y =  0.6,
    label = paste0("% with Y[i] > X[i]: ", round(100 * rando_frac, 1), "%"),
    
    hjust = 0,
    vjust = 1,
    size = 5
  )+ylim(0,0.6)+xlim(0,0.6)

grid.arrange(rando_plot_uns, immun_plot_uns, pleio_plot_uns, devel_plot_uns,nrow = 1)



######################## test of proportions ########################

K <- c(
  Developmental = nrow(above_line_devel),
  Immune = nrow(above_line_immun),
  Pleiotropic = nrow(above_line_pleio),
  Random = nrow(above_line_rando)
)

N = c(
  nrow(Developmental_uns),
  nrow(Immune_uns),
  nrow(Pleiotropic_uns),
  nrow(Random_uns)
)

prop.test(K, N)

pairwise.prop.test(K, N, p.adjust.method = "BH")



################################## mixed model and ICC ##################################
#mixed model
library(lme4)

Developmental_clean = Developmental_uns %>%
  filter(
    is.finite(W_domain),
    is.finite(W_gene),
    !is.na(protein)
  )


max(Developmental_clean$W_domain)
#997.9788

#remove large outliers
Developmental_clean = Developmental_clean[-which(Developmental_clean$W_domain> 2),]


dev_model  = lmer(W_domain ~ W_gene + (1 | protein), data = Developmental_clean)
imm_model  = lmer(W_domain ~ W_gene + (1 | protein), data = Immune_uns)
ple_model  = lmer(W_domain ~ W_gene + (1 | protein), data = Pleiotropic_uns)
ran_model  = lmer(W_domain ~ W_gene + (1 | protein), data = Random_uns)


############### check to see if the difference in fixed effect (W_gene) is sig ###########


all_dat = rbind(
  Developmental_clean,
  Immune_uns,
  Pleiotropic_uns,
  Random_uns
)
library(lme4)
library(lmerTest)

combined_model = lmerTest::lmer(
  W_domain ~ W_gene * group + (1 | protein),
  data = all_dat
)
summary(combined_model)
anova(combined_model)
library(emmeans )
emtrends(combined_model, pairwise ~ group, var = "W_gene")
detach("package:lmerTest", unload = TRUE)

######################## plot variance components Fig S5 ########################


get_icc <- function(model) {
  vc <- as.data.frame(VarCorr(model))
  
  var_protein  <- vc$vcov[vc$grp == "protein"]
  var_residual <- vc$vcov[vc$grp == "Residual"]
  
  icc <- var_protein / (var_protein + var_residual)
  
  data.frame(
    var_protein = var_protein,
    var_residual = var_residual,
    ICC = icc
  )
}


get_icc(ple_model)
get_icc(imm_model)
get_icc(dev_model)
get_icc(ran_model)

myresiduals = c(0.001735991, 0.001944873, 0.001073543, 0.001065273)
myprot = c(0.000441505, 0.001298271, 0.000666132, 0.004365288)
mynames = c('Pleiotropic', "Immunity", "Developmental", "Random")

res_prot = bind_rows(
  data.frame(
    group = mynames,
    variance = myresiduals,
    type = "Residual"
  ),
  data.frame(
    group = mynames,
    variance = myprot,
    type = "Protein"
  )
)

res_prot$group = factor(res_prot$group, levels = c('Pleiotropic', 'Immunity', 'Developmental', 'Random'))
ggplot(res_prot, aes(x = group, y = variance, color = type, group = type)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("blue", "red"),
    labels = c("Gene", "Residual")
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "Variance",
    color = ""
  ) +
  theme(
    axis.text.x = element_text(size = 13, angle = 90),
    axis.text.y = element_text(size = 13)
  )


########################## bootstraping ICC (Figure S4) ##########################

library(lme4)

boot_icc = function(model) {
  vc = as.data.frame(VarCorr(model))
  vp = vc$vcov[vc$grp == "protein"]
  vr = vc$vcov[vc$grp == "Residual"]
  vp / (vp + vr)
}

#do bootstaping using simulated values. The values are simulated using models  (e.g., pleio_kinase_model)
boot_pleio = bootMer(ple_model, boot_icc, nsim = 1000)
boot_immun = bootMer(imm_model, boot_icc, nsim = 1000)
boot_devel = bootMer(dev_model, boot_icc, nsim = 1000)
boot_rando = bootMer(ran_model, boot_icc, nsim = 1000)

diff_imm_pleio_icc = boot_immun$t - boot_pleio$t
diff_dev_pleio_icc = boot_devel$t - boot_pleio$t
diff_ran_pleio_icc = boot_rando$t - boot_pleio$t


mean(diff_imm_pleio_icc > 0)*100
mean(diff_dev_pleio_icc > 0)*100
mean(diff_ran_pleio_icc > 0)*100


library(tidyverse)

icc_diff_df = tibble(
  value = c(
    diff_imm_pleio_icc,
    diff_dev_pleio_icc,
    diff_ran_pleio_icc
  ),
  comparison = c(
    rep("Immune - Pleiotropic", length(diff_imm_pleio_icc)),
    rep("Developmental - Pleiotropic", length(diff_dev_pleio_icc)),
    rep("Random - Pleiotropic", length(diff_ran_pleio_icc))
  )
)



ggplot(icc_diff_df, aes(x = value,
                        fill = comparison,
                        color = comparison)) +
  geom_density(alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic(base_size = 13) +
  labs(
    x = "Bootstrap ICC Difference",
    y = "Density"
  ) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(
    values = c(
      "Immune - Pleiotropic" = "#65c2a5",
      "Developmental - Pleiotropic" = "#fc8d62",
      "Random - Pleiotropic" = "grey"
    )
  ) +
  scale_color_manual(
    values = c(
      "Immune - Pleiotropic" = "#65c2a5",
      "Developmental - Pleiotropic" = "#fc8d62",
      "Random - Pleiotropic" = "grey"
    )
  )

##############################################################################
