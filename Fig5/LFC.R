
Ortho_Imm_lfc = read.csv('Ortho_Imm_LFC.csv',row.names = 1)
Ortho_Imm_lfc = na.omit(Ortho_Imm_lfc)


Ortho_ple_lfc = read.csv('Ortho_Plei_LFC.csv',row.names = 1)
Ortho_ple_lfc = na.omit(Ortho_ple_lfc)



library(ggplot2)

# Calculate row means
imm = na.omit(as.numeric(rowMeans(Ortho_Imm_lfc[, 2:5], na.rm = TRUE)))
ple = na.omit(as.numeric(rowMeans(Ortho_ple_lfc[, 2:5], na.rm = TRUE)))

# Combine into a data frame
df = data.frame(
  value = c(imm, ple),
  group = c(rep("Immunity", length(imm)),
            rep("Pleiotropy", length(ple)))
)

# Plot
ggplot(df, aes(x = value, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  theme_classic() +
  labs(
    x = "Average log2FC across orthologs",
    y = "Density"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Immunity" = "#65c2a5",
      "Pleiotropy" = "#fc8d62"
    )
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "Immunity" = "#65c2a5",
      "Pleiotropy" = "#fc8d62"
    )
  )



Ortho_ran_lfc = read.csv('Ortho_random_LFC.csv',row.names = 1)

Ortho_ran_lfc = na.omit(Ortho_ran_lfc)


ran = na.omit(as.numeric(rowMeans(Ortho_ran_lfc[, 2:5], na.rm = TRUE)))


df = data.frame(
  value = c(imm, ple, ran),
  group = c(rep("Immunity", length(imm)),
            rep("Pleiotropy", length(ple)),
            rep("Random", length(ran)))
)


ggplot(df, aes(x = value, fill = group)) +
  geom_density(alpha = 1) +
  facet_wrap(~group, ncol = 3) +
  theme_classic(base_size = 14)+
  scale_color_manual(
    values = c(
      name = NULL,
      "Immunity" = "#65c2a5",
      "Pleiotropy" = "#fc8d62",
      "Random" = "grey50"
    )
  ) +
  scale_fill_manual(
    values = c(
      name = NULL,
      "Immunity" = "#65c2a5",
      "Pleiotropy" = "#fc8d62",
      "Random" = "grey50"
    )
  )


ggplot(df, aes(x = value, color = group)) +
  geom_density(linewidth = 1) +
  coord_cartesian(xlim = c(-2, 2)) +
  theme_classic(base_size = 14)+
  scale_color_manual(
    values = c(
      "Immunity" = "#65c2a5",
      "Pleiotropy" = "#fc8d62",
      "Random" = "grey50"
    ),
    labels = c(
      "Immunity" = "Non-plei\nImmune",
      "Pleiotropy" = "Pleiotropic",
      "Random" = "Random"
    )
  )  +
  scale_fill_manual(
    values = c(
      name = NULL,
      "Immunity" = "#65c2a5",
      "Pleiotropy" = "#fc8d62",
      "Random" = "grey50"
    )
  )+labs( x = "Average log2FC across orthologs", y = "Density" )+theme(legend.title = element_blank())

table(df$group)

