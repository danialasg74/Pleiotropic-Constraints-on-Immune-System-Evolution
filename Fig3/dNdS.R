################################### finding immunity classes ################################### 

library(jsonlite)
library(dplyr)
library(stringr)
library(ggplot2)


Sackton = readxl::read_excel('Sackton_2007.xls', sheet = 'Sheet1')
William = readxl::read_excel('Williams_2023.xlsx', sheet = 'Supplemental_Table_2')
Earlyyy = readxl::read_excel('Early_2017.xlsx', sheet = 'Sheet1')

Sackton = Sackton[,c(1,2,3,7)]
William = William[,c(1,4,5)]
Earlyyy = Earlyyy[,c(3,7,1)]

names(Sackton)[1] = 'gene'
names(William)[1] = 'gene'
names(Earlyyy)[1] = 'gene'


names(Sackton)[4] = 'class'
names(William)[3] = 'class'
names(Earlyyy)[2] = 'class'


Sackton = Sackton[,c(1,4,3)]
William = William[,c(1,3,2)]


names(Sackton)[3] = 'name'
names(William)[3] = 'name'
names(Earlyyy)[3] = 'name'

nrow(Sackton) + nrow(William) + nrow(Earlyyy)
df = rbind(Sackton, William, Earlyyy)

df = df[df$class != "NA", ]
df$class =  sub(" .*", "", df$class)
df$class = tolower(df$class)  
unique(df$class)

df = subset(df, !duplicated(df$gene))
names(df)


mart = read.delim('mart_export.txt')
mart = mart[,c(1,3)]
names(mart) = c('gene', 'ncbi')
mart = subset(mart, !duplicated(mart$gene))

dfncbi = merge(df,mart, by = 'gene' )

setdiff(df$gene,mart$gene )

ncbi = read.delim('ncbi_dataset.tsv')

ncbi = ncbi[,c(8,12,14)]
ncbi = ncbi[order(-ncbi$Protein.length), ]
ncbi = subset(ncbi, !duplicated(ncbi$Gene.ID))
ncbi = ncbi[,c(1,2)]

names(ncbi) = c('ncbi', 'protein')

dfncbi = merge(dfncbi, ncbi, by ='ncbi')
dfncbi = dfncbi[,c(5,3,4)]


ncbi_ann = read.delim('ncbi_droso_ann.tsv')
ncbi_ann = ncbi_ann[,c(6,7,12)]
names(ncbi_ann)[3] = 'protein'



################################### immunity dn/ds ################################### 


setwd('immun/')


# read JSON file
immun = fromJSON("dN_and_ds_plus_non_domains.txt")

# convert the nested list to a flat data frame
immun = do.call(rbind, lapply(names(immun), function(name) {
  cbind(file = name, as.data.frame(immun[[name]]))
}))
rownames(immun) <- NULL

immun$protein =str_extract(immun$file, "NP_[0-9]+\\.[0-9]+")

immun$domain <- ifelse(
  str_detect(immun$file, "/non_domain/"),
  "non_domain",
  str_extract(immun$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")
)


# for rows without a domain (just the main protein), set domain = "none"
immun$domain[is.na(immun$domain)] <- "none"

# reorder columns
immun = immun[, c("protein", "domain", "dN", "dS")]

immun = merge(immun, ncbi_ann, by = 'protein')
immun$W = (immun$dN)/(immun$dS)



##### make df of dS for domains genes and non-domains
immun = immun %>%
  mutate(region = case_when(
    domain == "non_domain" ~ "Non-domain",
    domain == "none"       ~ "Gene",
    !is.na(domain)         ~ "Domain",
    TRUE                   ~ NA_character_
  ))

immun$group = 'Immune'

immun_uns = subset(immun, immun$dS < 3)

immun_non_domain = immun_uns %>%
  filter(region == "Non-domain") %>%
  select(protein, W_non_domain = W)


immun_domain_non_domain = immun_uns %>%
  filter(region == "Domain") %>%
  left_join(immun_non_domain, by = "protein")




ggplot(immun_domain_non_domain, aes(x = W_non_domain, y = W)) +
  geom_point(size = 0.5, alpha = 1, color = "#65c2a5") +
  geom_smooth(method = "lm", color = "#65c2a5") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1.2) +
  annotate("text", x = 0.3, y = 0.45, label = "y = x", color = "black", size = 8, fontface = "italic") +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1), limits = c(0, 0.6)) +
  xlab("dN/dS Gene") +
  ylab("dN/dS Domain") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  )+ggtitle("Non-plei Immune")



IMM = immun_domain_non_domain[,c(1,2,7,10)]


################################### Pleio dn/ds ################################### 

setwd('/pleio/')


# read JSON file
pleio = fromJSON("dN_and_ds_plus_non_domains.txt")

# convert the nested list to a flat data frame
pleio = do.call(rbind, lapply(names(pleio), function(name) {
  cbind(file = name, as.data.frame(pleio[[name]]))
}))
rownames(pleio) <- NULL

pleio$protein =str_extract(pleio$file, "NP_[0-9]+\\.[0-9]+")

pleio$domain <- ifelse(
  str_detect(pleio$file, "/non_domain/"),
  "non_domain",
  str_extract(pleio$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")
)


# for rows without a domain (just the main protein), set domain = "none"
pleio$domain[is.na(pleio$domain)] <- "none"

# reorder columns
pleio = pleio[, c("protein", "domain", "dN", "dS")]

pleio = merge(pleio, ncbi_ann, by = 'protein')
pleio$W = (pleio$dN)/(pleio$dS)



##### make df of dS for domains genes and non-domains
pleio = pleio %>%
  mutate(region = case_when(
    domain == "non_domain" ~ "Non-domain",
    domain == "none"       ~ "Gene",
    !is.na(domain)         ~ "Domain",
    TRUE                   ~ NA_character_
  ))

pleio$group = 'pleioe'

pleio_uns = subset(pleio, pleio$dS < 3)

pleio_non_domain = pleio_uns %>%
  filter(region == "Non-domain") %>%
  select(protein, W_non_domain = W)


pleio_domain_non_domain = pleio_uns %>%
  filter(region == "Domain") %>%
  left_join(pleio_non_domain, by = "protein")




PLE = pleio_domain_non_domain[,c(1,2,7,10)]



IMM_mer  = merge(dfncbi, IMM, by = 'protein')
PLE_mer  = merge(dfncbi, PLE, by = 'protein')

IMM_mer$kind = 'Immunity'
PLE_mer$kind = 'Pleiotropic'


Final_dn_ds = rbind(IMM_mer,PLE_mer )
length(unique(Final_dn_ds$protein))


unique(Final_dn_ds$class)

Final_dn_ds$class = factor(Final_dn_ds$class, levels = c("recognition",
                                                         "signaling"  ,
                                                         "effector"   ,
                                                         'anti-viral' ,
                                                         'other'      ))

DF = Final_dn_ds

rm(list = setdiff(ls(), "DF"))




#################################### non-domain

DF_nondomain = subset(DF, !duplicated(DF$protein))

ggplot(DF_nondomain, aes( x = reorder(class, W_non_domain, FUN = function(x) median(x, na.rm = TRUE)), y = W_non_domain, fill = class ))+
  geom_violin(trim = FALSE, alpha = 1) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(
    values = c(
      "recognition" = "burlywood",
      "signaling"   = "palevioletred1",
      "effector"    = "lightblue",
      'anti-viral'  = 'turquoise',
      'other'       = 'beige'
    )
  ) +
  theme_classic()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11),
    axis.title.y = element_text(size = 13)
  ) +
  labs(title = "", x = "", y = "dN/dS (omega)")+
  labs(title = "", x = "", y = "dN/dS (omega)")+
  annotate(
    "text",
    x = 1,
    y = 2.3,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 2,
    y = 2.3,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 3,
    y = 2.3,
    label = "ab",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 4,
    y = 2.3,
    label = "b",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 5,
    y = 2.3,
    label = "ab",
    colour = "black",
    size = 5
  )


x = c(DF_nondomain$W_non_domain)
grp = DF_nondomain$class


kruskal.test(x ~ grp)


pairwise.wilcox.test(x, grp, p.adjust.method = "BH", exact = FALSE, na.action = na.omit)


aggregate(x ~ grp, FUN = median)

# Other a
# Signaling a
# Recognition ab
# Effector b
# Anti-viral ab

