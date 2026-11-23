
################################################  Import packages ######################################## 

########################################################################################################## 

library(jsonlite)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(emmeans)
library(patchwork)


mybasesize = 10

################################################  Pleiotropy ############################################## 

########################################################################################################### 



setwd('.../dnds/pleio')


# read JSON file
data <- fromJSON("dN_and_ds_plus_non_domains.txt")

# convert the nested list to a flat data frame
df <- do.call(rbind, lapply(names(data), function(name) {
  cbind(file = name, as.data.frame(data[[name]]))
}))
rownames(df) <- NULL

# extract protein name 
df$protein <- str_extract(df$file, "NP_[0-9]+\\.[0-9]+")

# getthe domains
df$domain <- ifelse(
  str_detect(df$file, "/non_domain/"),
  "non_domain",
  str_extract(df$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")
)


# for rows without a domain (just the main protein), set domain = "none"
df$domain[is.na(df$domain)] <- "none"

# reorder columns
df <- df[, c("protein", "domain", "dN", "dS")]

# import annotation
ncbi_ann = read.delim('/Users/danial/Documents/postdoc_research/Domains/ncbi_droso_ann.tsv')
ncbi_ann = ncbi_ann[,c(6,7,12)]
names(ncbi_ann)[3] = 'protein'

# dN/dS
df = merge(df, ncbi_ann, by = 'protein')
df$W = (df$dN)/(df$dS)


# unsaturated
df_uns = subset(df, df$dS<3)


# supplemental figure
pleio_d_g_n = ggplot(df_uns, aes(x = kind, y = W, fill = kind)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = c("Gene" = "#fc8d62", "Domain" = "red" , "Non-domain" = "black")) +
  labs(y = "dN/dS", x = "", title = "Pleiotropic") +
  theme_classic(base_size =mybasesize) +
  theme(legend.position = "none",plot.title = element_text(size = 13))



# Domain vs gene vs non-domain stat
x = c(df_uns$W)
grp = c(df_uns$kind)

kruskal.test(x ~ grp)


pairwise.wilcox.test(x, grp, p.adjust.method = "BH", exact = FALSE, na.action = na.omit)


#what proportion survives dS<3
nrow(df_uns)/nrow(df)


# plot domain versus non-domain

non_domain = df_uns %>%
  filter(kind == "Non-domain") %>%
  select(protein, W_non_domain = W)


domain_non_domain = df_uns %>%
  filter(kind == "Domain") %>%
  left_join(non_domain, by = "protein")


ggplot(domain_non_domain, aes(x = W_non_domain, y = W)) +
  geom_point() +
  theme_classic() +
  labs(
    x = expression(omega["Non-domain"]),
    y = expression(omega["Domain"])
  )

above_line_pleio = domain_non_domain[domain_non_domain$W > domain_non_domain$W_non_domain, ]

pleio_frac = nrow(above_line_pleio)/nrow(domain_non_domain)


mypleio = domain_non_domain

ple_model  = lmer(W ~ W_non_domain + (1 | protein), data = mypleio)

ple_int   = fixef(ple_model)[["(Intercept)"]]
ple_slope = fixef(ple_model)[["W_non_domain"]]


plot_pleio = ggplot(mypleio, aes(x = W_non_domain, y = W)) +
  geom_point(size = 0.5, alpha = 1, color = "#fc8d62") +
  geom_abline(
    intercept = ple_int,
    slope = ple_slope,
    color = "#fc8d62",
    linewidth = 1
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 2.2, by = 0.5), limits = c(0, 2.2)) +
  xlab("dN/dS non-domain") +
  ylab("dN/dS domain") +
  theme_classic(base_size = mybasesize) +
  ggtitle("Pleiotropic")#+ annotate( "text", x = 0.3, y = 0.7, label = paste0("% with Y[i] > X[i]: ", round(100 * pleio_frac, 1), "%"), hjust = 0, vjust = 1, size = 5 )





rm(list = setdiff(ls(), c("mypleio", "pleio_frac","plot_pleio", "pleio_d_g_n","mybasesize")))



################################################  Immunity ################################################ 

########################################################################################################### 


setwd('.../dnds/immun')


# read JSON file
data <- fromJSON("dN_and_ds_plus_non_domains.txt")

# convert the nested list to a flat data frame
df <- do.call(rbind, lapply(names(data), function(name) {
  cbind(file = name, as.data.frame(data[[name]]))
}))
rownames(df) <- NULL

# extract protein names
df$protein <- str_extract(df$file, "NP_[0-9]+\\.[0-9]+")

# get domains
df$domain <- ifelse(
  str_detect(df$file, "/non_domain/"),
  "non_domain",
  str_extract(df$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")
)



# for rows without a domain (just the main protein), set domain = "none"
df$domain[is.na(df$domain)] <- "none"

# reorder columns
df <- df[, c("protein", "domain", "dN", "dS")]

# annotation
ncbi_ann = read.delim('/Users/danial/Documents/postdoc_research/Domains/ncbi_droso_ann.tsv')
ncbi_ann = ncbi_ann[,c(6,7,12)]
names(ncbi_ann)[3] = 'protein'


df = merge(df, ncbi_ann, by = 'protein')
df$W = (df$dN)/(df$dS)


# unsaturated
df_uns = subset(df, df$dS<3)

# supplemental
immun_d_g_n = ggplot(df_uns, aes(x = kind, y = W, fill = kind)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = c("Gene" = "#65c2a5", "Domain" = "darkgreen" , "Non-domain" = "black")) +
  labs(y = "", x = "", title = "Non-plei \nImmune") +
  theme_classic(base_size =mybasesize) +
  theme(legend.position = "none",plot.title = element_text(size = 13))

# stat
x = c(df_uns$W)
grp = c(df_uns$kind)

kruskal.test(x ~ grp)


pairwise.wilcox.test(x, grp, p.adjust.method = "BH", exact = FALSE, na.action = na.omit)


#what proportion survives dS<3
nrow(df_uns)/nrow(df)



# domain versus non-domain:


non_domain = df_uns %>%
  filter(kind == "Non-domain") %>%
  select(protein, W_non_domain = W)


domain_non_domain = df_uns %>%
  filter(kind == "Domain") %>%
  left_join(non_domain, by = "protein")


ggplot(domain_non_domain, aes(x = W_non_domain, y = W)) +
  geom_point() +
  theme_classic() +
  labs(
    x = expression(omega["Non-domain"]),
    y = expression(omega["Domain"])
  )

above_line_immun = domain_non_domain[domain_non_domain$W > domain_non_domain$W_non_domain, ]

immun_frac = nrow(above_line_immun)/nrow(domain_non_domain)




myimmun = domain_non_domain

imm_model  = lmer(W ~ W_non_domain + (1 | protein), data = myimmun)

imm_int   = fixef(imm_model)[["(Intercept)"]]
imm_slope = fixef(imm_model)[["W_non_domain"]]

plot_immun = ggplot(domain_non_domain, aes(x = W_non_domain, y = W)) +
  geom_point(size = 0.5, alpha = 1, color = "#65c2a5") +
  geom_abline(
    intercept = imm_int,
    slope = imm_slope,
    color = "#65c2a5",
    linewidth = 1
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 2.2, by = 0.5), limits = c(0, 2.2)) +
  xlab("dN/dS non-domain") +
  ylab("") +
  theme_classic(base_size = mybasesize)+ggtitle("Non-plei \nImmune")#+ annotate( "text", x = 0.3, y = 0.7, label = paste0("% with Y[i] > X[i]: ", round(100 * immun_frac, 1), "%"), hjust = 0, vjust = 1, size = 5 )




rm(list = setdiff(ls(), c("mypleio", "pleio_frac", "plot_pleio",
                          "myimmun", "immun_frac", "plot_immun", "immun_d_g_n", "pleio_d_g_n","mybasesize")))



################################################  Developmental ################################################ 

################################################################################################################ 



setwd('.../dnds/devel')


# read JSON file
data <- fromJSON("dN_and_ds_plus_non_domains.txt")

# convert the nested list to a flat data frame
df <- do.call(rbind, lapply(names(data), function(name) {
  cbind(file = name, as.data.frame(data[[name]]))
}))
rownames(df) <- NULL



df$protein <- str_extract(df$file, "NP_[0-9]+\\.[0-9]+")


df$domain <- ifelse(
  str_detect(df$file, "/non_domain/"),
  "non_domain",
  str_extract(df$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")
)


# for rows without a domain (just the main protein), set domain = "none"
df$domain[is.na(df$domain)] <- "none"

# reorder columns
df <- df[, c("protein", "domain", "dN", "dS")]

ncbi_ann = read.delim('/Users/danial/Documents/postdoc_research/Domains/ncbi_droso_ann.tsv')
ncbi_ann = ncbi_ann[,c(6,7,12)]
names(ncbi_ann)[3] = 'protein'


df = merge(df, ncbi_ann, by = 'protein')
df$W = (df$dN)/(df$dS)

df_uns = subset(df, df$dS<3)


head(df_uns)

df_clean = df_uns %>%
  filter(
    is.finite(W),
    !is.na(protein)
  )


max(df_clean$W)
#1153
df_clean = df_clean[-which(df_clean$W> 2),]
df_uns = df_clean
max(df_uns$W)


devel_d_g_n = ggplot(df_uns, aes(x = kind, y = W, fill = kind)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = c("Gene" = "#8da0cb", "Domain" = "blue" , "Non-domain" = "black")) +
  labs(y = "", x = "", title = "Non-plei \nDevelopmental") +
  theme_classic(base_size =mybasesize) +
  theme(legend.position = "none",plot.title = element_text(size = 13))


x = c(df_uns$W)
grp = c(df_uns$kind)

kruskal.test(x ~ grp)


pairwise.wilcox.test(x, grp, p.adjust.method = "BH", exact = FALSE, na.action = na.omit)


#what proportion survives dS<3
nrow(df_uns)/nrow(df)






non_domain = df_uns %>%
  filter(kind == "Non-domain") %>%
  select(protein, W_non_domain = W)


domain_non_domain = df_uns %>%
  filter(kind == "Domain") %>%
  left_join(non_domain, by = "protein")


ggplot(domain_non_domain, aes(x = W_non_domain, y = W)) +
  geom_point() +
  theme_classic() +
  labs(
    x = expression(omega["Non-domain"]),
    y = expression(omega["Domain"])
  )


above_line_devel = domain_non_domain[domain_non_domain$W > domain_non_domain$W_non_domain, ]

devel_frac = nrow(above_line_devel)/nrow(domain_non_domain)




mydevel = domain_non_domain

dev_model  = lmer(W ~ W_non_domain + (1 | protein), data = mydevel)

dev_int   = fixef(dev_model)[["(Intercept)"]]
dev_slope = fixef(dev_model)[["W_non_domain"]]

plot_devel = ggplot(domain_non_domain, aes(x = W_non_domain, y = W)) +
  geom_point(size = 0.5, alpha = 1, color = "#8da0cb") +
  geom_abline(
    intercept = dev_int,
    slope = dev_slope,
    color = "#8da0cb",
    linewidth = 1
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 2.2, by = 0.5), limits = c(0, 2.2)) +
  xlab("dN/dS non-domain") +
  ylab("") +
  theme_classic(base_size = mybasesize)+ggtitle("Non-plei \nDevelopmental")#+annotate( "text", x = 0.3, y = 0.7, label = paste0("% with Y[i] > X[i]: ", round(100 * devel_frac, 1), "%"), hjust = 0, vjust = 1, size = 5 )




rm(list = setdiff(ls(), c("mypleio", "pleio_frac", "plot_pleio", 
                          "myimmun", "immun_frac","plot_immun",
                          "mydevel", "devel_frac","plot_devel",
                          "pleio_d_g_n","immun_d_g_n","devel_d_g_n","mybasesize")))






################################################  random ################################################ 

######################################################################################################### 




setwd('.../dnds/rando')


# read JSON file
data <- fromJSON("dN_and_ds_plus_non_domains.txt")

# convert the nested list to a flat data frame
df <- do.call(rbind, lapply(names(data), function(name) {
  cbind(file = name, as.data.frame(data[[name]]))
}))
rownames(df) <- NULL

# extract proteins
df$protein <- str_extract(df$file, "NP_[0-9]+\\.[0-9]+")

df$domain <- ifelse(
  str_detect(df$file, "/non_domain/"),
  "non_domain",
  str_extract(df$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")
)

# for rows without a domain (just the main protein), set domain = "none"
df$domain[is.na(df$domain)] <- "none"

# reorder columns
df <- df[, c("protein", "domain", "dN", "dS")]

ncbi_ann = read.delim('/Users/danial/Documents/postdoc_research/Domains/ncbi_droso_ann.tsv')
ncbi_ann = ncbi_ann[,c(6,7,12)]
names(ncbi_ann)[3] = 'protein'


df = merge(df, ncbi_ann, by = 'protein')
df$W = (df$dN)/(df$dS)

df_uns = subset(df, df$dS<3)


head(df_uns)

df_clean = df_uns %>%
  filter(
    is.finite(W),
    !is.na(protein)
  )


max(df_clean$W)
df_clean = df_clean[-which(df_clean$W> 3),]
df_uns = df_clean
max(df_uns$W)


rando_d_g_n = ggplot(df_uns, aes(x = kind, y = W, fill = kind)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = c("Gene" = "grey", "Domain" = "grey30" , "Non-domain" = "black")) +
  labs(y = "", x = "", title = "Random") +
  theme_classic(base_size =mybasesize) +
  theme(legend.position = "none",plot.title = element_text(size = 13))


x = c(df_uns$W)
grp = c(df_uns$kind)

kruskal.test(x ~ grp)


pairwise.wilcox.test(x, grp, p.adjust.method = "BH", exact = FALSE, na.action = na.omit)


#what proportion survives dS<3
nrow(df_uns)/nrow(df)

head(df_uns)




non_domain = df_uns %>%
  filter(kind == "Non-domain") %>%
  select(protein, W_non_domain = W)


domain_non_domain = df_uns %>%
  filter(kind == "Domain") %>%
  left_join(non_domain, by = "protein")


ggplot(domain_non_domain, aes(x = W_non_domain, y = W)) +
  geom_point() +
  theme_classic() +
  labs(
    x = expression(omega["Non-domain"]),
    y = expression(omega["Domain"])
  )


above_line_rando = domain_non_domain[domain_non_domain$W > domain_non_domain$W_non_domain, ]

rando_frac = nrow(above_line_rando)/nrow(domain_non_domain)


myrando = domain_non_domain

ran_model  = lmer(W ~ W_non_domain + (1 | protein), data = myrando)

ran_int   = fixef(ran_model)[["(Intercept)"]]
ran_slope = fixef(ran_model)[["W_non_domain"]]


plot_rando = ggplot(domain_non_domain, aes(x = W_non_domain, y = W)) +
  geom_point(size = 0.5, alpha = 1, color = "grey") +
  geom_abline(
    intercept = ran_int,
    slope = ran_slope,
    color = "grey",
    linewidth = 1
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 2.2, by = 0.5), limits = c(0, 2.2)) +
  xlab("dN/dS non-domain") +
  ylab("") +
  theme_classic(base_size = mybasesize)+ggtitle("Random") #+ annotate( "text", x = 0.3, y = 0.7, label = paste0("% with Y[i] > X[i]: ", round(100 * rando_frac, 1), "%"), hjust = 0, vjust = 1, size = 5 )


model_rando = lm(W ~ W_non_domain, data = domain_non_domain)
coef(model_rando)
summary(model_rando)$r.squared







rm(list = setdiff(ls(), c("mypleio", "pleio_frac", "plot_pleio", 
                          "myimmun", "immun_frac","plot_immun",
                          "mydevel", "devel_frac","plot_devel",
                          "myrando", "rando_frac","plot_rando",
                          "pleio_d_g_n","immun_d_g_n","devel_d_g_n","rando_d_g_n","mybasesize")))


grid.arrange(plot_pleio, plot_immun, plot_devel, plot_rando,nrow = 1)

grid.arrange(pleio_d_g_n, immun_d_g_n, devel_d_g_n, rando_d_g_n,nrow = 1)



##### Lines from linear mixed model


lines_lmxmodel <- data.frame(
  Group = c("Developmental", "Immunity", "Pleiotropy", "Random"),
  intercept = c(
    0.004095,
    0.02895,
    0.01746,
    0.01604
    
  ),
  slope = c(
    0.284053,
    0.30006,
    0.13738,
    0.32548
    
  )
)

lines_lmxmodel$Group = factor(lines_lmxmodel$Group, c('Random', 'Immunity', 'Developmental', 'Pleiotropy') )

lines_lmxmodel_plot = ggplot(lines_lmxmodel) +
  geom_abline(
    aes(
      intercept = intercept,
      slope = slope,
      color = Group
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = c(
      "Developmental" = "#8da0cb",
      "Immunity"      = "#65c2a5",
      "Pleiotropy"    = "#fc8d62",
      "Random"        = "grey"
    )
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_classic(base_size = 10) +
  labs(
    x = "dN/dS (non-domain)",
    y = "Predicted dN/dS (domain)",
    color = NULL
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1)
  )


lines_lmxmodel_plot_nolegend = ggplot(lines_lmxmodel) +
  geom_abline(
    aes(
      intercept = intercept,
      slope = slope,
      color = Group
    ),
    linewidth = 1
  ) +
  scale_color_manual(
    values = c(
      "Developmental" = "#8da0cb",
      "Immunity"      = "#65c2a5",
      "Pleiotropy"    = "#fc8d62",
      "Random"        = "grey"
    )
  ) +
  theme_classic(base_size = 10) +
  labs(
    x = "dN/dS (non-domain)",
    y = "Predicted dN/dS (domain)",
    color = NULL
  ) +theme(legend.position = "")+geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  scale_y_continuous( breaks = seq(0, 1, by = 0.2), limits = c(0, 1), expand = c(0, 0) ) + scale_x_continuous( breaks = seq(0, 2.2, by = 0.5), limits = c(0, 2.2), expand = c(0, 0) )

grid.arrange(plot_pleio, plot_immun, plot_devel, plot_rando,lines_lmxmodel_plot_nolegend,nrow = 1)



########################################


