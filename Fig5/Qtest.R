
#################### import files 

Ortho_Imm_LFC = read.csv("Ortho_Imm_LFC.csv",row.names = 1)
Ortho_Imm_LSE = read.csv("Ortho_Imm_LFCSE.csv",row.names = 1)





#H₀: LFC is homogeneous across the four species.
#H₁: LFC is heterogeneous 

species = c("mel", "sec", "sim", "yak")
# Keep only genes with LFC and SE values for all four species
complete_genes = complete.cases(
  Ortho_Imm_LFC[, species],
  Ortho_Imm_LSE[, species]
)

Imm_LFC = Ortho_Imm_LFC[complete_genes, ]
Imm_LSE = Ortho_Imm_LSE[complete_genes, ]

nrow(Ortho_Imm_LFC)
nrow(Imm_LFC)



Imm_LFC$Q = NA
Imm_LFC$df = NA
Imm_LFC$pvalue = NA
Imm_LFC$weighted_mean = NA

for (i in 1:nrow(Imm_LFC)) {
  
  lfc = as.numeric(
    Imm_LFC[i, c("mel", "sec", "sim", "yak")]
  )
  
  se = as.numeric(
    Imm_LSE[i, c("mel", "sec", "sim", "yak")]
  )
  
  # Weighted mean (more precise numbers get more weight)
  w = 1 / se^2
  
  weighted_mean = sum(w * lfc) / sum(w)
  Imm_LFC$weighted_mean[i] = weighted_mean
  # Q test: SUM [1/variance * (LFC - weighted mean LFC)^2]
  Q = sum(w * (lfc - weighted_mean)^2)
  
  # Degrees of freedom
  df = length(lfc) - 1
  
  # Heterogeneity p-value
  pvalue = pchisq(Q, df = df, lower.tail = FALSE)
  
  # Store results
  Imm_LFC$Q[i] = Q
  Imm_LFC$df[i] = df
  Imm_LFC$pvalue[i] = pvalue
}


Imm_LFC$FDR = p.adjust(Imm_LFC$pvalue, method = "BH")
Imm_LFC$heterogeneous = Imm_LFC$FDR < 0.05
table(Imm_LFC$heterogeneous)

rm(list = setdiff(ls(), "Imm_LFC"))

ggplot(Imm_LFC, aes(weighted_mean, Q)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_classic() +
  labs(
    x = "Weighted mean log2 fold change",
    y = "Cochran's Q"
  )

Imm_LFC$abs_weighted_mean = abs(Imm_LFC$weighted_mean)

ggplot(Imm_LFC,
       aes(abs_weighted_mean, Q)) +

  geom_point(alpha = 0.4, size = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  theme_classic() +
  labs(
    x = "Absolute weighted mean log2 fold change",
    y = "Cochran's Q"
  )

"Response magnitude was positively associated with cross-species heterogeneity."
cor.test(
  abs(Imm_LFC$weighted_mean),
  Imm_LFC$Q,
  method = "spearman"
)

#################### pleio ########################


Ortho_Ple_LFC = read.csv("Ortho_Plei_LFC.csv",row.names = 1)
Ortho_Ple_LSE = read.csv("Ortho_Plei_LFCSE.csv",row.names = 1)


species = c("mel", "sec", "sim", "yak")
# Keep only genes with LFC and SE values for all four species
complete_genes = complete.cases(
  Ortho_Ple_LFC[, species],
  Ortho_Ple_LSE[, species]
)

Ple_LFC = Ortho_Ple_LFC[complete_genes, ]
Ple_LSE = Ortho_Ple_LSE[complete_genes, ]

nrow(Ortho_Ple_LFC)
nrow(Ple_LFC)



Ple_LFC$Q = NA
Ple_LFC$df = NA
Ple_LFC$pvalue = NA


for (i in 1:nrow(Ple_LFC)) {
  
  lfc = as.numeric(
    Ple_LFC[i, c("mel", "sec", "sim", "yak")]
  )
  
  se = as.numeric(
    Ple_LSE[i, c("mel", "sec", "sim", "yak")]
  )
  
  # Weighted mean (more precise numbers get more weight)
  w = 1 / se^2
  
  weighted_mean = sum(w * lfc) / sum(w)

  
  # Q test: SUM [1/variance * (LFC - weighted mean LFC)^2]
  Q = sum(w * (lfc - weighted_mean)^2)
  
  # Degrees of freedom
  df = length(lfc) - 1
  
  # Heterogeneity p-value
  pvalue = pchisq(Q, df = df, lower.tail = FALSE)
  
  # Store results
  Ple_LFC$Q[i] = Q
  Ple_LFC$df[i] = df
  Ple_LFC$pvalue[i] = pvalue
}


Ple_LFC$FDR = p.adjust(Ple_LFC$pvalue, method = "BH")
Ple_LFC$heterogeneous = Ple_LFC$FDR < 0.05
table(Ple_LFC$heterogeneous)

rm(list = setdiff(ls(), c("Imm_LFC", "Ple_LFC")))

#################### random ########################

Ortho_Ran_LFC = read.csv("Ortho_random_LFC.csv",row.names = 1)
Ortho_Ran_LSE = read.csv("Ortho_random_LFCSE.csv",row.names = 1)


species = c("mel", "sec", "sim", "yak")
# Keep only genes with LFC and SE values for all four species
complete_genes = complete.cases(
  Ortho_Ran_LFC[, species],
  Ortho_Ran_LSE[, species]
)

Ran_LFC = Ortho_Ran_LFC[complete_genes, ]
Ran_LSE = Ortho_Ran_LSE[complete_genes, ]

nrow(Ortho_Ran_LFC)
nrow(Ran_LFC)



Ran_LFC$Q = NA
Ran_LFC$df = NA
Ran_LFC$pvalue = NA

for (i in 1:nrow(Ran_LFC)) {
  
  lfc = as.numeric(
    Ran_LFC[i, c("mel", "sec", "sim", "yak")]
  )
  
  se = as.numeric(
    Ran_LSE[i, c("mel", "sec", "sim", "yak")]
  )
  
  # Weighted mean (more precise numbers get more weight)
  w = 1 / se^2
  
  weighted_mean = sum(w * lfc) / sum(w)
  
  # Q test: SUM [1/variance * (LFC - weighted mean LFC)^2]
  Q = sum(w * (lfc - weighted_mean)^2)
  
  # Degrees of freedom
  df = length(lfc) - 1
  
  # Heterogeneity p-value
  pvalue = pchisq(Q, df = df, lower.tail = FALSE)
  
  # Store results
  Ran_LFC$Q[i] = Q
  Ran_LFC$df[i] = df
  Ran_LFC$pvalue[i] = pvalue
}


Ran_LFC$FDR = p.adjust(Ran_LFC$pvalue, method = "BH")
Ran_LFC$heterogeneous = Ran_LFC$FDR < 0.05
table(Ran_LFC$heterogeneous)


rm(list = setdiff(ls(), c("Imm_LFC", "Ple_LFC","Ran_LFC")))



table(Imm_LFC$heterogeneous)
table(Ple_LFC$heterogeneous)
table(Ran_LFC$heterogeneous)





heterogeneity_table = matrix(
  c(
    4, 326,
    3, 286,
    3, 887
  ),
  nrow = 3,
  byrow = TRUE
)

rownames(heterogeneity_table) = c(
  "Immunity",
  "Pleiotropic",
  "Random"
)

colnames(heterogeneity_table) = c(
  "Heterogeneous",
  "Not_heterogeneous"
)

heterogeneity_table
fisher.test(heterogeneity_table)



#######overall degree of LFC heterogeneity differ among immunity, pleiotropic, and random genes?#####

###plotting Q 


Imm_LFC$Group = "Immunity"
Ple_LFC$Group = "Pleiotropic"
Ran_LFC$Group = "Random"

Q_all = rbind(
  Imm_LFC[, c("Q", "Group")],
  Ple_LFC[, c("Q", "Group")],
  Ran_LFC[, c("Q", "Group")]
)


library(ggplot2)

ggplot(Q_all, aes(x = Q, fill = Group)) +
  geom_density(alpha = 0.4) +
  theme_classic() +
  xlab("Cochran's Q") +
  ylab("Density")

ggplot(Q_all, aes(x = Q, fill = Group)) +
  geom_density(alpha = 1) +
  scale_fill_manual(values = c(
    "Immunity" = "#65c2a5",
    "Pleiotropic" = "#fc8d62",
    "Random" = "grey"
  )) +
  theme_classic(base_size = 13) +
  xlab("Cochran's Q") +
  ylab("Density")+xlim(0,10)



ggplot(Q_all, aes(x = Q, fill = Group)) +
  geom_density(alpha = 0.6) +
  scale_x_log10() +
  theme_classic(base_size = 14) +
  xlab("Cochran's Q (log10 scale)") +
  ylab("Density")+
  scale_fill_manual(values = c(
    "Immunity" = "#65c2a5",
    "Pleiotropic" = "#fc8d62",
    "Random" = "grey"
  ))+theme(legend.title = element_blank())
summary(Imm_LFC$Q)
summary(Ple_LFC$Q)
summary(Ran_LFC$Q)

#Median Q: Immunity (1.34)> Random (0.97)> Pleiotropic (0.76)

kruskal.test(Q ~ Group, data = Q_all)
pairwise.wilcox.test(
  Q_all$Q,
  Q_all$Group,
  p.adjust.method = "BH"
)



medians = data.frame(
  Group = c("Immunity","Pleiotropic","Random"),
  Median = c(
    median(Imm_LFC$Q),
    median(Ple_LFC$Q),
    median(Ran_LFC$Q)
  )
)

ggplot(Q_all, aes(x = Q, fill = Group)) +
  geom_density(alpha = 0.6) +
  geom_vline(data = medians,
             aes(xintercept = Median, colour = Group),
             linetype = 2,
             linewidth = 0.8,
             show.legend = FALSE) +
  scale_x_log10() +
  theme_classic(base_size = 14) +
  xlab("Cochran's Q") +
  ylab("Density") +
  scale_fill_manual(values = c(
    "Immunity" = "#66c2a5",
    "Pleiotropic" = "#fc8d62",
    "Random" = "grey70"
  )) +
  scale_colour_manual(values = c(
    "Immunity" = "#66c2a5",
    "Pleiotropic" = "#fc8d62",
    "Random" = "grey40"
  )) +
  annotate("text",
           x = 50,
           y = 0.72,
           label = "P < 0.001",
           hjust = 1,
           size = 4.5) +
  theme(legend.title = element_blank())



library(ggridges)
ggplot(Q_all, aes(x = Q, y = Group, fill = Group)) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 1.2,
    color = "black",
    rel_min_height = 0.001
  ) +
  geom_vline(
    data = medians,
    aes(xintercept = Median, colour = Group),
    linetype = 2,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_x_log10() +
  scale_fill_manual(values = c(
    "Immunity" = "#66c2a5",
    "Random" = "grey70",
    "Pleiotropic" = "#fc8d62"
  )) +
  scale_colour_manual(values = c(
    "Immunity" = "#66c2a5",
    "Random" = "grey40",
    "Pleiotropic" = "#fc8d62"
  )) +
  scale_y_discrete(
    limits = c("Pleiotropic", "Random", "Immunity"),
    labels = c("Pleiotropic", "Random", "Non-plei\nImmune"),
    expand = expansion(add = c(0.5, 2.5))
  )  +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 17) +
  xlab("Cochran's Q (log10 scale)") +
  ylab(NULL)  +
  theme(legend.position = "none")+
  annotate(
    "label",
    x = 0.001,
    y = 4.8,
    hjust = 0,
    size = 3.8,
    family = "mono",
    label =
      "Pairwise Wilcoxon
Non-plei-Imm vs Plei   <0.001
Non-plei-Imm vs Random   <0.001
Plei vs Random    0.027"
  )

  