
################################### R/C ################################### 


setwd('/Users/danial/Documents/postdoc_research/Domains/My_approach/seq_3_groups/main/updated/AllDomains/sixspeices/Radical_over_Cons')
pleio = read.csv('pleio_all.csv',row.names = 1)
immun = read.csv('immun_all.csv',row.names = 1)
devel = read.csv('devel_all.csv',row.names = 1)



################################### finding immunity classes ################################### 


#we used previous literature (Early et al. 2017; Sackton et al. 2007) 
#and supplemental Table 2 of Williams et al. (2023) 
#to categorize pleiotropic and non-pleiotropic immune genes into 
#recognition, signaling, effectors, antiviral, and other immune genes. 

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


mart = read.delim('../mart_export.txt')
mart = mart[,c(1,3)]
names(mart) = c('gene', 'ncbi')
mart = subset(mart, !duplicated(mart$gene))

dfncbi = merge(df,mart, by = 'gene' )

setdiff(df$gene,mart$gene )

ncbi = read.delim('../ncbi_dataset.tsv')

ncbi = ncbi[,c(8,12,14)]
ncbi = ncbi[order(-ncbi$Protein.length), ]
ncbi = subset(ncbi, !duplicated(ncbi$Gene.ID))
ncbi = ncbi[,c(1,2)]

names(ncbi) = c('ncbi', 'protein')

dfncbi = merge(dfncbi, ncbi, by ='ncbi')
dfncbi = dfncbi[,c(5,3,4)]


rm(list = setdiff(ls(), c("dfncbi","pleio",'immun', 'devel')))






pleio_uniq = subset(pleio, !duplicated(pleio$protein))
immun_uniq = subset(immun, !duplicated(immun$protein))
devel_uniq = subset(devel, !duplicated(devel$protein))
pleio_uniq = pleio_uniq[,c(1,9)]
immun_uniq = immun_uniq[,c(1,9)]
devel_uniq = devel_uniq[,c(1,9)]
pleio_uniq$kind = 'pleio'
immun_uniq$kind = 'immun'
devel_uniq$kind = 'devel'
#Adding Devel makes no difference
GENES = rbind(pleio_uniq, immun_uniq )
GENES_class = merge(GENES, dfncbi, by = 'protein')


GENES_class$class = factor(GENES_class$class, levels = c("recognition",
                                                         "signaling"  ,
                                                         "effector"   ,
                                                         'anti-viral' ,
                                                         'other'      ))


ggplot(GENES_class, aes(x = class, y = RdC_pro, fill = class)) +
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
    y = 2.7,
    label = "ab",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 2,
    y = 2.7,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 3,
    y = 2.7,
    label = "b",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 4,
    y = 2.7,
    label = "ab",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 5,
    y = 2.7,
    label = "ab",
    colour = "black",
    size = 5
  )
  
  
  
#recognition  "ab"
#signaling    "a"
#effector     "b"
#anti-viral   "ab"
#other        "ab"
  
  
  
x = c(GENES_class$RdC_pro)
grp = GENES_class$class


kruskal.test(x ~ grp)


pairwise.wilcox.test(x, grp, p.adjust.method = "BH", exact = FALSE, na.action = na.omit)

aggregate(x ~ grp, FUN = median)



subpleio = pleio[,c(1,2,9,10)]
subimmun = immun[,c(1,2,9,10)]
subdevel = devel[,c(1,2,9,10)]
subpleio$kind = 'pleio'
subimmun$kind = 'immun'

DOMS = rbind(subpleio, subimmun )
DOMS_class = merge(DOMS, dfncbi, by = 'protein')

DOMS_class$class = factor(DOMS_class$class, levels = c("recognition",
                                                         "signaling"  ,
                                                         "effector"   ,
                                                         'anti-viral' ,
                                                         'other'      ))

DOMS_class$RdC_dom
ggplot(DOMS_class, aes(x = RdC_pro, y = RdC_dom)) +
  geom_point(aes(color = class), size = 1, alpha = 0.5) +
  #geom_smooth(method = "lm", color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "black", size = 1.2) +
  annotate("text", x =1.5, y = 2, label = "y = x",
           color = "black", size = 8, fontface = "italic") +
  #scale_y_continuous(breaks = seq(0, 0.6, by = 0.1),
  #                   limits = c(0, 0.6)) +
  xlab("R/C Protein") +
  ylab("R/C Domain") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  ggtitle("")+scale_color_manual(values = c(      "recognition" = "burlywood",
                                                  "signaling"   = "palevioletred1",
                                                  "effector"    = "lightblue",
                                                  'anti-viral'  = 'turquoise',
                                                  'other'       = 'beige'
  ))
above_line = DOMS_class[DOMS_class$RdC_dom > DOMS_class$RdC_pro, ]
prop_df =  data.frame(
  class = names(table(DOMS_class$class)),
  proportion = as.numeric(table(above_line$class)[names(table(DOMS_class$class))]) /
    as.numeric(table(DOMS_class$class))
)


library(ggplot2)
prop_df = prop_df %>%
  mutate(class = reorder(class, -proportion))
ggplot(prop_df, aes(x = class, y = proportion)) +
  geom_col(fill = "steelblue") +
  theme_classic(base_size = 15) +
  ylab("Proportion above y = x") +
  xlab("") +
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(prop_df, aes(x = class, y = proportion, fill = class)) +
  geom_col(color = "black") +
  theme_classic(base_size = 15) +
  ylab("Proportion above y = x") +
  xlab("") +
  ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "recognition" = "burlywood",
    "signaling"   = "palevioletred1",
    "effector"    = "lightblue",
    "anti-viral"  = "turquoise",
    "other"       = "beige"
  ))+theme(legend.position  = "")








################ domain ################
ggplot(DOMS_class, aes(x = class, y = RdC_pro, fill = class)) +
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
    y = 2.7,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 2,
    y = 2.7,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 3,
    y = 2.7,
    label = "b",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 4,
    y = 2.7,
    label = "b",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 5,
    y = 2.7,
    label = "b",
    colour = "black",
    size = 5
  )
#recognition  "a"
#signaling    "a"
#effector     "b"
#anti-viral   "b"
#other        "b"

x = c(DOMS_class$RdC_pro)
grp = DOMS_class$class


kruskal.test(x ~ grp)

 
pairwise.wilcox.test(x, grp, p.adjust.method = "BH", exact = FALSE, na.action = na.omit)



aggregate(x ~ grp, FUN = median)
