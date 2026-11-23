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

################################### immunity R/C ################################### 

setwd('immun')


x = readLines("RoverC.txt")

# Remove null lines
x = x[x != "null/null/null/null"]

# Split on /
x = strsplit(x, "/")

# Protein name is always first
protein = sapply(x, `[`, 1)

# Region: protein-level rows have only 4 fields
region = sapply(x, function(z) {
  if (length(z) == 4) "none" else z[2]
})

# R, C, Ratio depend on whether there are 4 or 5 fields
R = as.numeric(sapply(x, function(z) {
  if (length(z) == 4) z[2] else z[3]
}))

C = as.numeric(sapply(x, function(z) {
  if (length(z) == 4) z[3] else z[4]
}))

Ratio = as.numeric(sapply(x, function(z) {
  if (length(z) == 4) z[4] else z[5]
}))

df = data.frame(
  protein = protein,
  domain = region,
  R = R,
  C = C,
  Ratio = Ratio
)


#### remove those with dS<3
mydnds = read.csv("immun_dnds.csv",row.names = 1)
#all bellow 3?
sum(mydnds$dS < 3) == nrow(mydnds)
mydnds$w = mydnds$dN/ mydnds$dS
mydnds = mydnds[,c("protein","domain","w")]

df = df %>%
  left_join(
    mydnds %>% select(protein, domain, w),
    by = c("protein", "domain")
  ) %>%
  filter(!is.na(w))

df = subset(df, select = -w)
###########################

mynondomain = df %>%
  filter(domain == "non_domain") %>%
  select(protein, Ratio_non_domain = Ratio)


domain_non_domain = df %>%
  filter(
    !domain %in% c("none", "non_domain", "non_domain_all")
  ) %>%
  left_join(mynondomain, by = "protein")

IMM = domain_non_domain[,c(1,2,5,6)]

rm(list = setdiff(ls(), c("IMM", "dfncbi")))



################################### Pleiotropy R/C ################################### 

setwd('pleio')


x = readLines("RoverC.txt")

# Remove null lines
x = x[x != "null/null/null/null"]

# Split on /
x = strsplit(x, "/")

# Protein name is always first
protein = sapply(x, `[`, 1)

# Region: protein-level rows have only 4 fields
region = sapply(x, function(z) {
  if (length(z) == 4) "none" else z[2]
})

# R, C, Ratio depend on whether there are 4 or 5 fields
R = as.numeric(sapply(x, function(z) {
  if (length(z) == 4) z[2] else z[3]
}))

C = as.numeric(sapply(x, function(z) {
  if (length(z) == 4) z[3] else z[4]
}))

Ratio = as.numeric(sapply(x, function(z) {
  if (length(z) == 4) z[4] else z[5]
}))

df = data.frame(
  protein = protein,
  domain = region,
  R = R,
  C = C,
  Ratio = Ratio
)


#### remove those with dS<3
mydnds = read.csv("pleio_dnds.csv",row.names = 1)
#all bellow 3?
sum(mydnds$dS < 3) == nrow(mydnds)
mydnds$w = mydnds$dN/ mydnds$dS
mydnds = mydnds[,c("protein","domain","w")]

df = df %>%
  left_join(
    mydnds %>% select(protein, domain, w),
    by = c("protein", "domain")
  ) %>%
  filter(!is.na(w))

df = subset(df, select = -w)
###########################

mynondomain = df %>%
  filter(domain == "non_domain") %>%
  select(protein, Ratio_non_domain = Ratio)


domain_non_domain = df %>%
  filter(
    !domain %in% c("none", "non_domain", "non_domain_all")
  ) %>%
  left_join(mynondomain, by = "protein")

PLE = domain_non_domain[,c(1,2,5,6)]


rm(list = setdiff(ls(), c("IMM", "dfncbi", "PLE")))



IMM_mer  = merge(dfncbi, IMM, by = 'protein')
PLE_mer  = merge(dfncbi, PLE, by = 'protein')

IMM_mer$kind = 'Immunity'
PLE_mer$kind = 'Pleiotropic'


Final_dn_ds = rbind(IMM_mer,PLE_mer )

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

DF_nondomain$class
ggplot(DF_nondomain, aes( x = reorder(class, Ratio_non_domain, FUN = function(x) median(x, na.rm = TRUE)), y = Ratio_non_domain, fill = class ))+
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
  labs(title = "", x = "", y = "R/C")+
  annotate(
    "text",
    x = 1,
    y = 4,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 2,
    y = 4,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 3,
    y = 4,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 4,
    y = 4,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 5,
    y = 4,
    label = "a",
    colour = "black",
    size = 5
  )


x = c(DF_nondomain$Ratio_non_domain)
grp = DF_nondomain$class


kruskal.test(x ~ grp)


pairwise.wilcox.test(x, grp, p.adjust.method = "BH", exact = FALSE, na.action = na.omit)


aggregate(x ~ grp, FUN = median)


