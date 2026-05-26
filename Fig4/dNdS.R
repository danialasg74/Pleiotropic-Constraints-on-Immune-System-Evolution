library(jsonlite)
library(dplyr)
library(stringr)
library(ggplot2)


################################### finding immunity classes ################################### 

#we used previous literature (Early et al. 2017; Sackton et al. 2007) 
#and supplemental Table 2 of Williams et al. (2023) 
#to categorize pleiotropic and non-pleiotropic immune genes into 
#recognition, signaling, effectors, antiviral, and other immune genes. 

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






################################### immunity dn/ds ################################### 




ncbi_ann = read.delim('../ncbi_droso_ann.tsv')
ncbi_ann = ncbi_ann[,c(6,7,12)]
names(ncbi_ann)[3] = 'protein'
setwd('../immun')


# read JSON file
immun = fromJSON("dN_and_ds.txt")

# convert the nested list to a flat data frame
immun = do.call(rbind, lapply(names(immun), function(name) {
  cbind(file = name, as.data.frame(immun[[name]]))
}))
rownames(immun) = NULL

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

immun$group = 'Immune'

immun_uns = subset(immun, immun$dS < 3)

immun_uns$W_gene = NA

for (i in seq_len(nrow(immun_uns))) {
  PROT = subset(immun_uns, immun_uns$protein == immun_uns$protein[i])
  
  w_val <- subset(PROT, PROT$domain == "none")$W
  
  if (length(w_val) > 0) {
    # if there are multiple, just take the first; adjust if you want something else
    immun_uns$W_gene[immun_uns$protein == immun_uns$protein[i]] <- w_val[1]
  } else {
    cat("Skipped index:", i, " protein:", immun_uns$protein[i], "\n")
    next
  }
}

rm(PROT, i , w_val)


immun_uns_Dom = subset(immun_uns, immun_uns$kind=='Domain')


ggplot(immun_uns_Dom, aes(x = W_gene, y = W)) +
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
  )+ggtitle("Immunity")



IMM = immun_uns_Dom[,c(1,2,7,10)]


################################### Pleio dn/ds ################################### 

setwd('../pleio')


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

pleio$group = 'Pleio'

pleio_uns = subset(pleio, pleio$dS < 3)

pleio_uns$W_gene = NA

for (i in seq_len(nrow(pleio_uns))) {
  PROT = subset(pleio_uns, pleio_uns$protein == pleio_uns$protein[i])
  
  w_val <- subset(PROT, PROT$domain == "none")$W
  
  if (length(w_val) > 0) {
    # if there are multiple, just take the first; adjust if you want something else
    pleio_uns$W_gene[pleio_uns$protein == pleio_uns$protein[i]] <- w_val[1]
  } else {
    cat("Skipped index:", i, " protein:", pleio_uns$protein[i], "\n")
    next
  }
}

rm(PROT, i , w_val)

pleio_uns_Dom = subset(pleio_uns, pleio_uns$kind=='Domain')


PLE = pleio_uns_Dom[,c(1,2,7,10)]



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
ggplot(Final_dn_ds, aes(x = W_gene, y = W)) +
  geom_point(aes(color = class), size = 1, alpha = 0.5) +
  #geom_smooth(method = "lm", color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "black", size = 1.2) +
  annotate("text", x = 0.3, y = 0.45, label = "y = x",
           color = "black", size = 8, fontface = "italic") +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1),
                     limits = c(0, 0.6)) +
  xlab("dN/dS Gene") +
  ylab("dN/dS Domain") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  ggtitle("")+scale_color_manual(values = c(             "recognition" = "burlywood",
                                                         "signaling"   = "palevioletred1",
                                                         "effector"    = "lightblue",
                                                         'anti-viral'  = 'turquoise',
                                                         'other'       = 'beige'
                                                    ))
above_line = Final_dn_ds[Final_dn_ds$W > Final_dn_ds$W_gene, ]



prop_df =  data.frame(
  class = names(table(Final_dn_ds$class)),
  proportion = as.numeric(table(above_line$class)[names(table(Final_dn_ds$class))]) /
    as.numeric(table(Final_dn_ds$class))
)


library(ggplot2)

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

GENES = subset(Final_dn_ds, !duplicated(Final_dn_ds$protein))
GENES = GENES[,c(1,2,3,6,7)]
unique(GENES$class)

ggplot(GENES, aes(x = class, y = W_gene, fill = class)) +
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
  annotate(
    "text",
    x = 1,
    y = 0.88,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 2,
    y = 0.88,
    label = "b",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 3,
    y = 0.88,
    label = "a",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 4,
    y = 0.88,
    label = "ab",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 5,
    y = 0.88,
    label = "b",
    colour = "black",
    size = 5
  )

x = c(GENES$W_gene)
grp = GENES$class


kruskal.test(x ~ grp)


pairwise.wilcox.test(x, grp, p.adjust.method = "BH", exact = FALSE, na.action = na.omit)


#recognition  "a"
#signaling    "b"
#effector     "a"
#anti-viral   "ab"
#other        "b"

aggregate(x ~ grp, FUN = median)








############# domains #############

ggplot(Final_dn_ds, aes(x = class, y = W, fill = class)) +
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
  annotate(
    "text",
    x = 1,
    y = 0.88,
    label = "ab",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 2,
    y = 0.88,
    label = "d",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 3,
    y = 0.88,
    label = "ac",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 4,
    y = 0.88,
    label = "c",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 5,
    y = 0.88,
    label = "b",
    colour = "black",
    size = 5
  )



#recognition  "ab"
#signaling    "d"
#effector     "ac"
#anti-viral   "c"
#other        "b"



x = c(Final_dn_ds$W)
grp = Final_dn_ds$class


kruskal.test(x ~ grp)


pairwise.wilcox.test(x, grp, p.adjust.method = "BH", exact = FALSE, na.action = na.omit)


aggregate(x ~ grp, FUN = median)
