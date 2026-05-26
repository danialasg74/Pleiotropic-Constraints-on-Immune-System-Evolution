library(jsonlite)
library(dplyr)
library(stringr)
library(ggplot2)


################################### Processing R/C for the 4 classes of genes ###################################

###################################################################################################################

##### file for anntation ####

ann = read.delim('/Users/danial/Documents/postdoc_research/Domains/ncbi_droso_ann.tsv')
ann = ann[,c(6,8,12)]
names(ann)[3] = "protein"



####################### Processing Radical changes #######################



#################### pleiotropic genes ####################

setwd('../pleio')


# ==== CONFIG ====
path <- "protein_domain_ratios_all.txt"   # <-- set your file path

# ==== READ & PARSE ====
# Read raw lines to avoid accidental column splitting by spaces
lines <- read_lines(path)
# Drop blank lines if any
lines <- lines[nzchar(lines)]

# Split each line by "/" into parts
parts <- strsplit(lines, "/", fixed = TRUE)

# Extract fields:
# - protein: always the first part
# - ratio: always the last part
# - domain: everything between (usually length 0 or 1); NA if none
protein <- vapply(parts, function(x) x[1], character(1))
ratio   <- vapply(parts, function(x) x[length(x)], character(1))
domain  <- vapply(parts, function(x) {
  if (length(x) <= 2) return(NA_character_)
  paste(x[2:(length(x)-1)], collapse = "/")
}, character(1))

df <- tibble(
  protein = protein,
  domain  = domain,
  ratio   = suppressWarnings(as.numeric(ratio))
) %>%
  mutate(
    type = if_else(is.na(domain), "protein", "domain")
  ) %>%
  group_by(protein) %>%
  mutate(
    # take the ratio from the protein row (type == "protein")
    protein_ratio = first(ratio[is.na(domain)])
  ) %>%
  ungroup()


df= subset(df , df$type == 'domain')
names(df)[3] = 'domain_ratio'


Pleio = ggplot(df, aes(x = protein_ratio, y = domain_ratio)) +
  geom_point(alpha = 1,size = 0.1) +                        # scatter points
  theme_classic(base_size = 14) +                                # clean theme
  labs(
    x = "Protein",
    y = "Domain",
    title = "Proportion of Radical Changes"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )+geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.5)+xlim(0,1)+ylim(0,1)

subset(df, df$domain_ratio>0.75)



############## add pleiotropic dN/dS for filtering out genes with dS <3 ##############


# read JSON file
dnds <- fromJSON("dN_and_ds.txt")

# convert the nested list to a flat data frame
dnds <- do.call(rbind, lapply(names(dnds), function(name) {
  cbind(file = name, as.data.frame(dnds[[name]]))
}))
rownames(dnds) <- NULL

# extract protein and domain names using regex
library(stringr)

dnds$protein <- str_extract(dnds$file, "NP_[0-9]+\\.[0-9]+")
dnds$domain  <- str_extract(dnds$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")

# for rows without a domain (just the main protein), set domain = "none"
dnds$domain[is.na(dnds$domain)] <- "none"

# reorder columns
dnds <- dnds[, c("protein", "domain", "dN", "dS")]



dnds$kind  = ifelse(dnds$domain == "none", "Gene", "Domain")

dnds = subset(dnds, dnds$kind =='Domain')

dnds = dnds[,-ncol(dnds)]
df = df[,-4]
df_dnds = merge(df, dnds)
df_dnds_uns =  subset(df_dnds, df_dnds$dS < 3 )

Pleio_rad = df_dnds_uns

rm(list = setdiff(ls(), c("Pleio_rad")))



#################### Immune genes ####################


setwd('../immun')


# ==== CONFIG ====
path <- "protein_domain_ratios_all.txt"   # <-- set your file path

# ==== READ & PARSE ====
# Read raw lines to avoid accidental column splitting by spaces
lines <- read_lines(path)
# Drop blank lines if any
lines <- lines[nzchar(lines)]

# Split each line by "/" into parts
parts <- strsplit(lines, "/", fixed = TRUE)

# Extract fields:
# - protein: always the first part
# - ratio: always the last part
# - domain: everything between (usually length 0 or 1); NA if none
protein <- vapply(parts, function(x) x[1], character(1))
ratio   <- vapply(parts, function(x) x[length(x)], character(1))
domain  <- vapply(parts, function(x) {
  if (length(x) <= 2) return(NA_character_)
  paste(x[2:(length(x)-1)], collapse = "/")
}, character(1))

df <- tibble(
  protein = protein,
  domain  = domain,
  ratio   = suppressWarnings(as.numeric(ratio))
) %>%
  mutate(
    type = if_else(is.na(domain), "protein", "domain")
  ) %>%
  group_by(protein) %>%
  mutate(
    # take the ratio from the protein row (type == "protein")
    protein_ratio = first(ratio[is.na(domain)])
  ) %>%
  ungroup()


df= subset(df , df$type == 'domain')
names(df)[3] = 'domain_ratio'


Immun = ggplot(df, aes(x = protein_ratio, y = domain_ratio)) +
  geom_point(alpha = 1,size = 0.1) +                        # scatter points
  theme_classic(base_size = 14) +                                # clean theme
  labs(
    x = "Protein",
    y = "Domain",
    title = "Proportion of Radical Changes"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )+geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.5)+xlim(0,1)+ylim(0,1)

subset(df, df$domain_ratio>0.75)



############## add immune dN/dS for filtering out genes with dS <3 ##############

# read JSON file
dnds <- fromJSON("dN_and_ds.txt")

# convert the nested list to a flat data frame
dnds <- do.call(rbind, lapply(names(dnds), function(name) {
  cbind(file = name, as.data.frame(dnds[[name]]))
}))
rownames(dnds) <- NULL

# extract protein and domain names using regex
library(stringr)

dnds$protein <- str_extract(dnds$file, "NP_[0-9]+\\.[0-9]+")
dnds$domain  <- str_extract(dnds$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")

# for rows without a domain (just the main protein), set domain = "none"
dnds$domain[is.na(dnds$domain)] <- "none"

# reorder columns
dnds <- dnds[, c("protein", "domain", "dN", "dS")]



dnds$kind  = ifelse(dnds$domain == "none", "Gene", "Domain")

dnds = subset(dnds, dnds$kind =='Domain')

dnds = dnds[,-ncol(dnds)]
df = df[,-4]
df_dnds = merge(df, dnds)
df_dnds_uns =  subset(df_dnds, df_dnds$dS < 3 )

Immun_rad = df_dnds_uns


rm(list = setdiff(ls(), c("Pleio_rad","Immun_rad")))






#################### developmental genes ####################



setwd('../devel')


# ==== CONFIG ====
path <- "protein_domain_ratios_all.txt"   # <-- set your file path

# ==== READ & PARSE ====
# Read raw lines to avoid accidental column splitting by spaces
lines <- read_lines(path)
# Drop blank lines if any
lines <- lines[nzchar(lines)]

# Split each line by "/" into parts
parts <- strsplit(lines, "/", fixed = TRUE)

# Extract fields:
# - protein: always the first part
# - ratio: always the last part
# - domain: everything between (usually length 0 or 1); NA if none
protein <- vapply(parts, function(x) x[1], character(1))
ratio   <- vapply(parts, function(x) x[length(x)], character(1))
domain  <- vapply(parts, function(x) {
  if (length(x) <= 2) return(NA_character_)
  paste(x[2:(length(x)-1)], collapse = "/")
}, character(1))

df <- tibble(
  protein = protein,
  domain  = domain,
  ratio   = suppressWarnings(as.numeric(ratio))
) %>%
  mutate(
    type = if_else(is.na(domain), "protein", "domain")
  ) %>%
  group_by(protein) %>%
  mutate(
    # take the ratio from the protein row (type == "protein")
    protein_ratio = first(ratio[is.na(domain)])
  ) %>%
  ungroup()


df= subset(df , df$type == 'domain')
names(df)[3] = 'domain_ratio'


Devel = ggplot(df, aes(x = protein_ratio, y = domain_ratio)) +
  geom_point(alpha = 1,size = 0.1) +                        # scatter points
  theme_classic(base_size = 14) +                                # clean theme
  labs(
    x = "Protein",
    y = "Domain",
    title = "Proportion of Radical Changes"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )+geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.5)+xlim(0,1)+ylim(0,1)

subset(df, df$domain_ratio>0.75)





############## add developmental dN/dS for filtering out genes with dS <3 ##############


# read JSON file
dnds <- fromJSON("dN_and_ds.txt")

# convert the nested list to a flat data frame
dnds <- do.call(rbind, lapply(names(dnds), function(name) {
  cbind(file = name, as.data.frame(dnds[[name]]))
}))
rownames(dnds) <- NULL

# extract protein and domain names using regex
library(stringr)

dnds$protein <- str_extract(dnds$file, "NP_[0-9]+\\.[0-9]+")
dnds$domain  <- str_extract(dnds$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")

# for rows without a domain (just the main protein), set domain = "none"
dnds$domain[is.na(dnds$domain)] <- "none"

# reorder columns
dnds <- dnds[, c("protein", "domain", "dN", "dS")]



dnds$kind  = ifelse(dnds$domain == "none", "Gene", "Domain")

dnds = subset(dnds, dnds$kind =='Domain')

dnds = dnds[,-ncol(dnds)]
df = df[,-4]
df_dnds = merge(df, dnds)
df_dnds_uns =  subset(df_dnds, df_dnds$dS < 3 )

Devel_rad = df_dnds_uns


rm(list = setdiff(ls(), c("Pleio_rad","Immun_rad",'Devel_rad')))


############## random

setwd('../Random')

# ==== CONFIG ====
path <- "protein_domain_ratios_all.txt"   # <-- set your file path

# ==== READ & PARSE ====
# Read raw lines to avoid accidental column splitting by spaces
lines <- read_lines(path)
# Drop blank lines if any
lines <- lines[nzchar(lines)]

# Split each line by "/" into parts
parts <- strsplit(lines, "/", fixed = TRUE)

# Extract fields:
# - protein: always the first part
# - ratio: always the last part
# - domain: everything between (usually length 0 or 1); NA if none
protein <- vapply(parts, function(x) x[1], character(1))
ratio   <- vapply(parts, function(x) x[length(x)], character(1))
domain  <- vapply(parts, function(x) {
  if (length(x) <= 2) return(NA_character_)
  paste(x[2:(length(x)-1)], collapse = "/")
}, character(1))

df <- tibble(
  protein = protein,
  domain  = domain,
  ratio   = suppressWarnings(as.numeric(ratio))
) %>%
  mutate(
    type = if_else(is.na(domain), "protein", "domain")
  ) %>%
  group_by(protein) %>%
  mutate(
    # take the ratio from the protein row (type == "protein")
    protein_ratio = first(ratio[is.na(domain)])
  ) %>%
  ungroup()


df= subset(df , df$type == 'domain')
names(df)[3] = 'domain_ratio'


Devel = ggplot(df, aes(x = protein_ratio, y = domain_ratio)) +
  geom_point(alpha = 1,size = 0.1) +                        # scatter points
  theme_classic(base_size = 14) +                                # clean theme
  labs(
    x = "Protein",
    y = "Domain",
    title = "Proportion of Radical Changes"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )+geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.5)+xlim(0,1)+ylim(0,1)

subset(df, df$domain_ratio>0.75)




############## add random dN/dS for filtering out genes with dS <3 ##############

# read JSON file
dnds <- fromJSON("dN_and_ds.txt")

# convert the nested list to a flat data frame
dnds <- do.call(rbind, lapply(names(dnds), function(name) {
  cbind(file = name, as.data.frame(dnds[[name]]))
}))
rownames(dnds) <- NULL

# extract protein and domain names using regex
library(stringr)

dnds$protein <- str_extract(dnds$file, "NP_[0-9]+\\.[0-9]+")
dnds$domain  <- str_extract(dnds$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")

# for rows without a domain (just the main protein), set domain = "none"
dnds$domain[is.na(dnds$domain)] <- "none"

# reorder columns
dnds <- dnds[, c("protein", "domain", "dN", "dS")]



dnds$kind  = ifelse(dnds$domain == "none", "Gene", "Domain")

dnds = subset(dnds, dnds$kind =='Domain')

dnds = dnds[,-ncol(dnds)]
df = df[,-4]
df_dnds = merge(df, dnds)
df_dnds_uns =  subset(df_dnds, df_dnds$dS < 3 )

Rando_rad = df_dnds_uns


rm(list = setdiff(ls(), c("Pleio_rad","Immun_rad",'Devel_rad', 'Rando_rad')))



############# end of making dataframes for radical changes #############




####################### Processing Conserved changes #######################



#################### pleiotropic genes ####################


setwd('../pleio')

# ==== CONFIG ====
path <- "protein_domain_ratios_conserved_all.txt"   # <-- set your file path

# ==== READ & PARSE ====
# Read raw lines to avoid accidental column splitting by spaces
lines <- read_lines(path)
# Drop blank lines if any
lines <- lines[nzchar(lines)]

# Split each line by "/" into parts
parts <- strsplit(lines, "/", fixed = TRUE)

# Extract fields:
# - protein: always the first part
# - ratio: always the last part
# - domain: everything between (usually length 0 or 1); NA if none
protein <- vapply(parts, function(x) x[1], character(1))
ratio   <- vapply(parts, function(x) x[length(x)], character(1))
domain  <- vapply(parts, function(x) {
  if (length(x) <= 2) return(NA_character_)
  paste(x[2:(length(x)-1)], collapse = "/")
}, character(1))

df <- tibble(
  protein = protein,
  domain  = domain,
  ratio   = suppressWarnings(as.numeric(ratio))
) %>%
  mutate(
    type = if_else(is.na(domain), "protein", "domain")
  ) %>%
  group_by(protein) %>%
  mutate(
    # take the ratio from the protein row (type == "protein")
    protein_ratio = first(ratio[is.na(domain)])
  ) %>%
  ungroup()


df= subset(df , df$type == 'domain')
names(df)[3] = 'domain_ratio'


Pleio = ggplot(df, aes(x = protein_ratio, y = domain_ratio)) +
  geom_point(alpha = 1,size = 0.1) +                        # scatter points
  theme_classic(base_size = 14) +                                # clean theme
  labs(
    x = "Protein",
    y = "Domain",
    title = "Proportion of Radical Changes"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )+geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.5)+xlim(0,1)+ylim(0,1)

subset(df, df$domain_ratio>0.75)




############## add pleiotropic dN/dS for filtering out genes with dS <3 ##############


# read JSON file
dnds <- fromJSON("dN_and_ds.txt")

# convert the nested list to a flat data frame
dnds <- do.call(rbind, lapply(names(dnds), function(name) {
  cbind(file = name, as.data.frame(dnds[[name]]))
}))
rownames(dnds) <- NULL

# extract protein and domain names using regex
library(stringr)

dnds$protein <- str_extract(dnds$file, "NP_[0-9]+\\.[0-9]+")
dnds$domain  <- str_extract(dnds$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")

# for rows without a domain (just the main protein), set domain = "none"
dnds$domain[is.na(dnds$domain)] <- "none"

# reorder columns
dnds <- dnds[, c("protein", "domain", "dN", "dS")]



dnds$kind  = ifelse(dnds$domain == "none", "Gene", "Domain")

dnds = subset(dnds, dnds$kind =='Domain')

dnds = dnds[,-ncol(dnds)]
df = df[,-4]
df_dnds = merge(df, dnds)
df_dnds_uns =  subset(df_dnds, df_dnds$dS < 3 )

Pleio_cons = df_dnds_uns


rm(list = setdiff(ls(), c("Pleio_rad","Immun_rad",'Devel_rad', 'Rando_rad','Pleio_cons')))







#################### immune genes ####################

setwd('../immun')


# ==== CONFIG ====
path <- "protein_domain_ratios_conserved_all.txt"   # <-- set your file path

# ==== READ & PARSE ====
# Read raw lines to avoid accidental column splitting by spaces
lines <- read_lines(path)
# Drop blank lines if any
lines <- lines[nzchar(lines)]

# Split each line by "/" into parts
parts <- strsplit(lines, "/", fixed = TRUE)

# Extract fields:
# - protein: always the first part
# - ratio: always the last part
# - domain: everything between (usually length 0 or 1); NA if none
protein <- vapply(parts, function(x) x[1], character(1))
ratio   <- vapply(parts, function(x) x[length(x)], character(1))
domain  <- vapply(parts, function(x) {
  if (length(x) <= 2) return(NA_character_)
  paste(x[2:(length(x)-1)], collapse = "/")
}, character(1))

df <- tibble(
  protein = protein,
  domain  = domain,
  ratio   = suppressWarnings(as.numeric(ratio))
) %>%
  mutate(
    type = if_else(is.na(domain), "protein", "domain")
  ) %>%
  group_by(protein) %>%
  mutate(
    # take the ratio from the protein row (type == "protein")
    protein_ratio = first(ratio[is.na(domain)])
  ) %>%
  ungroup()


df= subset(df , df$type == 'domain')
names(df)[3] = 'domain_ratio'


Immun = ggplot(df, aes(x = protein_ratio, y = domain_ratio)) +
  geom_point(alpha = 1,size = 0.1) +                        # scatter points
  theme_classic(base_size = 14) +                                # clean theme
  labs(
    x = "Protein",
    y = "Domain",
    title = "Proportion of Radical Changes"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )+geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.5)+xlim(0,1)+ylim(0,1)

subset(df, df$domain_ratio>0.75)




############## add immune dN/dS for filtering out genes with dS <3 ##############

# read JSON file
dnds <- fromJSON("dN_and_ds.txt")

# convert the nested list to a flat data frame
dnds <- do.call(rbind, lapply(names(dnds), function(name) {
  cbind(file = name, as.data.frame(dnds[[name]]))
}))
rownames(dnds) <- NULL

# extract protein and domain names using regex
library(stringr)

dnds$protein <- str_extract(dnds$file, "NP_[0-9]+\\.[0-9]+")
dnds$domain  <- str_extract(dnds$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")

# for rows without a domain (just the main protein), set domain = "none"
dnds$domain[is.na(dnds$domain)] <- "none"

# reorder columns
dnds <- dnds[, c("protein", "domain", "dN", "dS")]



dnds$kind  = ifelse(dnds$domain == "none", "Gene", "Domain")

dnds = subset(dnds, dnds$kind =='Domain')

dnds = dnds[,-ncol(dnds)]
df = df[,-4]
df_dnds = merge(df, dnds)
df_dnds_uns =  subset(df_dnds, df_dnds$dS < 3 )

Immun_cons = df_dnds_uns
rm(list = setdiff(ls(), c("Pleio_rad","Immun_rad",'Devel_rad', 'Rando_rad','Pleio_cons', 'Immun_cons')))






#################### developmental genes ####################


setwd('../devel')


# ==== CONFIG ====
path <- "protein_domain_ratios_conserved_all.txt"   # <-- set your file path

# ==== READ & PARSE ====
# Read raw lines to avoid accidental column splitting by spaces
lines <- read_lines(path)
# Drop blank lines if any
lines <- lines[nzchar(lines)]

# Split each line by "/" into parts
parts <- strsplit(lines, "/", fixed = TRUE)

# Extract fields:
# - protein: always the first part
# - ratio: always the last part
# - domain: everything between (usually length 0 or 1); NA if none
protein <- vapply(parts, function(x) x[1], character(1))
ratio   <- vapply(parts, function(x) x[length(x)], character(1))
domain  <- vapply(parts, function(x) {
  if (length(x) <= 2) return(NA_character_)
  paste(x[2:(length(x)-1)], collapse = "/")
}, character(1))

df <- tibble(
  protein = protein,
  domain  = domain,
  ratio   = suppressWarnings(as.numeric(ratio))
) %>%
  mutate(
    type = if_else(is.na(domain), "protein", "domain")
  ) %>%
  group_by(protein) %>%
  mutate(
    # take the ratio from the protein row (type == "protein")
    protein_ratio = first(ratio[is.na(domain)])
  ) %>%
  ungroup()


df= subset(df , df$type == 'domain')
names(df)[3] = 'domain_ratio'


Devel = ggplot(df, aes(x = protein_ratio, y = domain_ratio)) +
  geom_point(alpha = 1,size = 0.1) +                        # scatter points
  theme_classic(base_size = 14) +                                # clean theme
  labs(
    x = "Protein",
    y = "Domain",
    title = "Proportion of Radical Changes"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )+geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.5)+xlim(0,1)+ylim(0,1)

subset(df, df$domain_ratio>0.75)





############## add developmental dN/dS for filtering out genes with dS <3 ##############


# read JSON file
dnds <- fromJSON("dN_and_ds.txt")

# convert the nested list to a flat data frame
dnds <- do.call(rbind, lapply(names(dnds), function(name) {
  cbind(file = name, as.data.frame(dnds[[name]]))
}))
rownames(dnds) <- NULL

# extract protein and domain names using regex
library(stringr)

dnds$protein <- str_extract(dnds$file, "NP_[0-9]+\\.[0-9]+")
dnds$domain  <- str_extract(dnds$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")

# for rows without a domain (just the main protein), set domain = "none"
dnds$domain[is.na(dnds$domain)] <- "none"

# reorder columns
dnds <- dnds[, c("protein", "domain", "dN", "dS")]



dnds$kind  = ifelse(dnds$domain == "none", "Gene", "Domain")

dnds = subset(dnds, dnds$kind =='Domain')

dnds = dnds[,-ncol(dnds)]
df = df[,-4]
df_dnds = merge(df, dnds)
df_dnds_uns =  subset(df_dnds, df_dnds$dS < 3 )
Devel_cons = df_dnds_uns

rm(list = setdiff(ls(), c("Pleio_rad","Immun_rad",'Devel_rad', 'Rando_rad', 'Pleio_cons', 'Immun_cons', 'Devel_cons')))



#################### random genes ####################

setwd('../Random')


# ==== CONFIG ====
path <- "protein_domain_ratios_conserved_all.txt"   # <-- set your file path

# ==== READ & PARSE ====
# Read raw lines to avoid accidental column splitting by spaces
lines <- read_lines(path)
# Drop blank lines if any
lines <- lines[nzchar(lines)]

# Split each line by "/" into parts
parts <- strsplit(lines, "/", fixed = TRUE)

# Extract fields:
# - protein: always the first part
# - ratio: always the last part
# - domain: everything between (usually length 0 or 1); NA if none
protein <- vapply(parts, function(x) x[1], character(1))
ratio   <- vapply(parts, function(x) x[length(x)], character(1))
domain  <- vapply(parts, function(x) {
  if (length(x) <= 2) return(NA_character_)
  paste(x[2:(length(x)-1)], collapse = "/")
}, character(1))

df <- tibble(
  protein = protein,
  domain  = domain,
  ratio   = suppressWarnings(as.numeric(ratio))
) %>%
  mutate(
    type = if_else(is.na(domain), "protein", "domain")
  ) %>%
  group_by(protein) %>%
  mutate(
    # take the ratio from the protein row (type == "protein")
    protein_ratio = first(ratio[is.na(domain)])
  ) %>%
  ungroup()


df= subset(df , df$type == 'domain')
names(df)[3] = 'domain_ratio'


Devel = ggplot(df, aes(x = protein_ratio, y = domain_ratio)) +
  geom_point(alpha = 1,size = 0.1) +                        # scatter points
  theme_classic(base_size = 14) +                                # clean theme
  labs(
    x = "Protein",
    y = "Domain",
    title = "Proportion of Radical Changes"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )+geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.5)+xlim(0,1)+ylim(0,1)

subset(df, df$domain_ratio>0.75)




############## add random dN/dS for filtering out genes with dS <3 ##############

# read JSON file
dnds <- fromJSON("dN_and_ds.txt")

# convert the nested list to a flat data frame
dnds <- do.call(rbind, lapply(names(dnds), function(name) {
  cbind(file = name, as.data.frame(dnds[[name]]))
}))
rownames(dnds) <- NULL

# extract protein and domain names using regex
library(stringr)

dnds$protein <- str_extract(dnds$file, "NP_[0-9]+\\.[0-9]+")
dnds$domain  <- str_extract(dnds$file, "PF[0-9]+\\.[0-9]+(_[0-9_]+)?")

# for rows without a domain (just the main protein), set domain = "none"
dnds$domain[is.na(dnds$domain)] <- "none"

# reorder columns
dnds <- dnds[, c("protein", "domain", "dN", "dS")]



dnds$kind  = ifelse(dnds$domain == "none", "Gene", "Domain")

dnds = subset(dnds, dnds$kind =='Domain')

dnds = dnds[,-ncol(dnds)]
df = df[,-4]
df_dnds = merge(df, dnds)
df_dnds_uns =  subset(df_dnds, df_dnds$dS < 3 )
Rando_cons = df_dnds_uns

rm(list = setdiff(ls(), c("Pleio_rad","Immun_rad",'Devel_rad', 'Rando_rad', 'Pleio_cons', 'Immun_cons', 'Devel_cons', 'Rando_cons')))


############ end of making dataframes for radical and conserved changes ##################





############ Plot R/C and calculate % with Y[i] > X[i] (Figure 1B) ##################



names(Pleio_rad)[3] = paste0(names(Pleio_rad)[3], '_rad')
names(Immun_rad)[3] = paste0(names(Immun_rad)[3], '_rad')
names(Devel_rad)[3] = paste0(names(Devel_rad)[3], '_rad')
names(Rando_rad)[3] = paste0(names(Rando_rad)[3], '_rad')
names(Pleio_rad)[4] = paste0(names(Pleio_rad)[4], '_rad')
names(Immun_rad)[4] = paste0(names(Immun_rad)[4], '_rad')
names(Devel_rad)[4] = paste0(names(Devel_rad)[4], '_rad')
names(Rando_rad)[4] = paste0(names(Rando_rad)[4], '_rad')

names(Pleio_cons)[3] = paste0(names(Pleio_cons)[3], '_cons')
names(Immun_cons)[3] = paste0(names(Immun_cons)[3], '_cons')
names(Devel_cons)[3] = paste0(names(Devel_cons)[3], '_cons')
names(Rando_cons)[3] = paste0(names(Rando_cons)[3], '_cons')
names(Pleio_cons)[4] = paste0(names(Pleio_cons)[4], '_cons')
names(Immun_cons)[4] = paste0(names(Immun_cons)[4], '_cons')
names(Devel_cons)[4] = paste0(names(Devel_cons)[4], '_cons')
names(Rando_cons)[4] = paste0(names(Rando_cons)[4], '_cons')

Pleio = merge(Pleio_rad, Pleio_cons, by = c('protein', 'domain'))
Immun = merge(Immun_rad, Immun_cons, by = c('protein', 'domain'))
Devel = merge(Devel_rad, Devel_cons, by = c('protein', 'domain'))
Rando = merge(Rando_rad, Rando_cons, by = c('protein', 'domain'))


Pleio = Pleio[,-c(5,6)]
Immun = Immun[,-c(5,6)]
Devel = Devel[,-c(5,6)]
Rando = Rando[,-c(5,6)]

names(Pleio)[7] = 'dN'
names(Immun)[7] = 'dN'
names(Devel)[7] = 'dN'
names(Rando)[7] = 'dN'

names(Pleio)[8] = 'dS'
names(Immun)[8] = 'dS'
names(Devel)[8] = 'dS'
names(Rando)[8] = 'dS'






Pleio$RdC_pro <- ifelse(Pleio$protein_ratio_cons == 0,0,Pleio$protein_ratio_rad / Pleio$protein_ratio_cons)
Immun$RdC_pro <- ifelse(Immun$protein_ratio_cons == 0,0,Immun$protein_ratio_rad / Immun$protein_ratio_cons)
Devel$RdC_pro <- ifelse(Devel$protein_ratio_cons == 0,0,Devel$protein_ratio_rad / Devel$protein_ratio_cons)
Rando$RdC_pro <- ifelse(Rando$protein_ratio_cons == 0,0,Rando$protein_ratio_rad / Rando$protein_ratio_cons)

Pleio$RdC_dom <- ifelse(Pleio$domain_ratio_cons == 0,0,Pleio$domain_ratio_rad / Pleio$domain_ratio_cons)
Immun$RdC_dom <- ifelse(Immun$domain_ratio_cons == 0,0,Immun$domain_ratio_rad / Immun$domain_ratio_cons)
Devel$RdC_dom <- ifelse(Devel$domain_ratio_cons == 0,0,Devel$domain_ratio_rad / Devel$domain_ratio_cons)
Rando$RdC_dom <- ifelse(Rando$domain_ratio_cons == 0,0,Rando$domain_ratio_rad / Rando$domain_ratio_cons)




Pleio = na.omit(Pleio)
above_line_pleio = Pleio[Pleio$RdC_dom > Pleio$RdC_pro, ]
Pleio_plot = ggplot(Pleio, aes(x = RdC_pro, y = RdC_dom)) +
  geom_point(alpha = 1,size = 0.1, color = "#fc8d62") + 
  geom_smooth(method = "lm", color = "#fc8d62")+
  annotate("text", x = 1, y = 4, label = paste0("% with Y[i] > X[i]: ", round(100 * nrow(above_line_pleio) / nrow(Pleio), 1), "%"), color = "black", size = 5)+
  theme_classic() +                                
  labs(
    x = "R/C Protein",
    y = "",
    title = "Pleiotropic"
  ) +
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    title = element_text(size = 14)
  )+geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7)+ylim(0,6)



Immun = na.omit(Immun)
above_line_immun = Immun[Immun$RdC_dom > Immun$RdC_pro, ]
Immun_plot = ggplot(Immun, aes(x = RdC_pro, y = RdC_dom)) +
  geom_point(alpha = 1,size = 0.1, color = "#65c2a5") +  
  annotate("text", x = 1, y = 4, label = paste0("% with Y[i] > X[i]: ", round(100 * nrow(above_line_immun) / nrow(Immun), 1), "%"), color = "black", size = 5)+
  geom_smooth(method = "lm", color = "#65c2a5")+
  theme_classic() +                                
  labs(
    x = "R/C Protein",
    y = "",
    title = "Immunity"
  ) +
  theme(  legend.title = element_blank(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          title = element_text(size = 14)
  )+geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7)+ylim(0,6)




Devel = na.omit(Devel)
above_line_devel = Devel[Devel$RdC_dom > Devel$RdC_pro, ]
Devel_plot = ggplot(Devel, aes(x = RdC_pro, y = RdC_dom)) +
  geom_point(alpha = 1,size = 0.1, color = "#8da0cb") +  
  geom_smooth(method = "lm", color = "#8da0cb")+
  annotate("text", x = 2, y = 4, label = paste0("% with Y[i] > X[i]: ", round(100 * nrow(above_line_devel) / nrow(Devel), 1), "%"), color = "black", size = 5)+
  theme_classic() +                                
  labs(
    x = "R/C Protein",
    y = "",
    title = "Developmental"
  ) +
  theme(    legend.title = element_blank(),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 14),
            legend.text = element_text(size = 14),
            title = element_text(size = 14)
            
  )+geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7)+ylim(0,6)




Rando = na.omit(Rando)
above_line_rando = Rando[Rando$RdC_dom > Rando$RdC_pro, ]
Rando_plot = ggplot(Rando, aes(x = RdC_pro, y = RdC_dom)) +
  geom_point(alpha = 1,size = 0.1, color = "grey") +  
  geom_smooth(method = "lm", color = "grey")+
  annotate("text", x = 1.2, y = 4, label = paste0("% with Y[i] > X[i]: ", round(100 * nrow(above_line_rando) / nrow(Rando), 1), "%"), color = "black", size = 5)+
  theme_classic() +                                
  labs(
    x = "R/C Protein",
    y = "R/C Domain",
    title = "Random"
  ) +
  theme(    legend.title = element_blank(),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 14),
            legend.text = element_text(size = 14),
            title = element_text(size = 14)
            
  )+geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7)+ylim(0,6)


library(gridExtra)
grid.arrange(Rando_plot, Immun_plot, Pleio_plot,Devel_plot,nrow  = 1)


######################## test of proportions ########################

K <- c(
  Developmental = nrow(above_line_devel),
  Immune = nrow(above_line_immun),
  Pleiotropic = nrow(above_line_pleio),
  Random = nrow(above_line_rando)
)

N = c(
  nrow(Devel),
  nrow(Immun),
  nrow(Pleio),
  nrow(Rando)
)

prop.test(K, N)

pairwise.prop.test(K, N, p.adjust.method = "BH")


#Developmental	a
#Pleiotropic	 ab
#Random	bc
#Immune	c



################################## mixed model and ICC ##################################
#mixed model
library(lme4)


dev_model  = lmer(RdC_dom ~ RdC_pro + (1 | protein), data = Devel)
imm_model  = lmer(RdC_dom ~ RdC_pro + (1 | protein), data = Immun)
ple_model  = lmer(RdC_dom ~ RdC_pro + (1 | protein), data = Pleio)
ran_model  = lmer(RdC_dom ~ RdC_pro + (1 | protein), data = Rando)

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


######################## plot variance components Fig S6 ########################


myresiduals = c(0.2269852, 0.1711354, 0.3590702, 0.1234386)
myprot = c(0.01637668, 0.02581867, 0.03652206, 0.01138143)
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
  theme_classic() +
  labs(
    x = "",
    y = "Variance",
    color = ""
  )+theme(axis.text.x = element_text(size = 13, angle =90),
          axis.text.y = element_text(size = 13))


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





##############################################################################
