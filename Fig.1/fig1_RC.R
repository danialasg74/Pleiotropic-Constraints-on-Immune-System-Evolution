library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(lmerTest)
library(patchwork)
library(emmeans)

######################### 
# Pleio
#########################

setwd('.../RC/pleio')

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


non_domains = df %>%
  filter(domain %in% c("non_domain", "non_domain_all")) %>%
  select(protein, domain, Ratio) %>%
  pivot_wider(
    names_from = domain,
    values_from = Ratio
  )

# calculate deviation from y=x
non_domains= non_domains %>%
  mutate(
    deviation = non_domain - non_domain_all
  )

# Deviation histogram
deviation_hist = ggplot(non_domains, aes(x = deviation)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    bins = 100
  ) +
  # geom_vline(
  #   xintercept = 0,
  #   linetype = "dashed",
  #   linewidth = 0.3
  # ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1)
  ) +
  theme_classic(base_size = 7) +
  labs(
    x = "Deviation from Y = X",
    y = "Proportion"
  )




Pleio_nondom = ggplot(non_domains, aes(x = non_domain_all, y = non_domain)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  coord_equal() +
  theme_classic(base_size = 13) +
  labs(
    x = "R/C ratio (non-domain-all)",
    y = "R/C ratio (non-domain)"
  )+ggtitle('Pleiotropy')

Pleio_nondom_inset = Pleio_nondom +
  inset_element(
    deviation_hist,
    left = 0.58,
    bottom = 0.55,
    right = 0.98,
    top = 0.98
  )



#### remove those with dS<3
mydnds = read.csv("/Users/danial/Documents/postdoc_research/Domains/My_approach/seq_3_groups/main/updated/AllDomains/sixspeices/a_a_changes/pleio/pleio_dnds.csv",row.names = 1)
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

#write
#write.csv(domain_non_domain, "/Users/danial/Documents/postdoc_research/Domains/My_approach/seq_3_groups/main/updated/AllDomains/sixspeices/a_a_changes/pleio/pleio_RC.csv")


# above y=x

pabove = nrow(domain_non_domain[domain_non_domain$Ratio > domain_non_domain$Ratio_non_domain, ])
ptotal = nrow(domain_non_domain)
above_line = pabove/ptotal * 100

# bindary assignment for each domain
binary = domain_non_domain %>%
  mutate(
    above = as.integer(Ratio > Ratio_non_domain)
  )

# Keep protein ID and binary assignment
binary = binary[, c("protein", "above")]
mypleio_binary = binary

# Remove domains where the comparison could not be made
mypleio_binary = mypleio_binary %>%
  filter(!is.na(above))

# name
mypleio_binary$group = "Pleiotropic"


## regression 
mypleio = domain_non_domain
ple_model  = lmer(Ratio ~ Ratio_non_domain + (1 | protein), data = domain_non_domain)

ple_int   = fixef(ple_model)[["(Intercept)"]]
ple_slope = fixef(ple_model)[["Ratio_non_domain"]]

Pleio_main = ggplot(
  domain_non_domain,
  aes(x = Ratio_non_domain, y = Ratio)
) +
  geom_abline(
    intercept = ple_int,
    slope = ple_slope,
    color = "#fc8d62",
    linewidth = 1
  ) +
  geom_point(size = 0.5, alpha = 1, color = "#fc8d62") +
  #geom_smooth(method = "lm", color = "#fc8d62") +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  theme_classic(base_size = 10) +
  labs(
    x = "R/C ratio (non-domain)",
    y = "R/C ratio (domain)"
  )+ coord_cartesian(xlim = c(0, 10), ylim = c(0, 10))#+ annotate( "text", x = 2.5, y = 7.5, label = paste0("%", round(above_line, 2)) )

rm(list = setdiff(ls(), c("Pleio_nondom", "Pleio_nondom_inset","Pleio_main","pabove","ptotal","mypleio",
                          "mypleio_binary")))






#############################
# immunity
##############################


setwd("..")
setwd("immun")




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




non_domains = df %>%
  filter(domain %in% c("non_domain", "non_domain_all")) %>%
  select(protein, domain, Ratio) %>%
  pivot_wider(
    names_from = domain,
    values_from = Ratio
  )


# calculate deviation from y=x
non_domains= non_domains %>%
  mutate(
    deviation = non_domain - non_domain_all
  )

# Deviation histogram
deviation_hist = ggplot(non_domains, aes(x = deviation)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    bins = 100
  ) +
  # geom_vline(
  #   xintercept = 0,
  #   linetype = "dashed",
  #   linewidth = 0.3
  # ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1)
  ) +
  theme_classic(base_size = 7) +
  labs(
    x = "Deviation from Y = X",
    y = "Proportion"
  )


immun_nondom = ggplot(non_domains, aes(x = non_domain_all, y = non_domain)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  coord_equal() +
  theme_classic(base_size = 13) +
  labs(
    x = "R/C ratio (non-domain-all)",
    y = ""
  )+ggtitle('Immunity')

immun_nondom_inset = immun_nondom +
  inset_element(
    deviation_hist,
    left = 0.58,
    bottom = 0.55,
    right = 0.98,
    top = 0.98
  )


#### remove those with dS<3
mydnds = read.csv("/Users/danial/Documents/postdoc_research/Domains/My_approach/seq_3_groups/main/updated/AllDomains/sixspeices/a_a_changes/immun/immun_dnds.csv",row.names = 1)
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



#write
#write.csv(domain_non_domain, "/Users/danial/Documents/postdoc_research/Domains/My_approach/seq_3_groups/main/updated/AllDomains/sixspeices/a_a_changes/immun/immun_RC.csv")



# above y=x
iabove = nrow(domain_non_domain[domain_non_domain$Ratio > domain_non_domain$Ratio_non_domain, ])
itotal = nrow(domain_non_domain)
above_line = iabove/itotal* 100

# bindary assignment for each domain
binary = domain_non_domain %>%
  mutate(
    above = as.integer(Ratio > Ratio_non_domain)
  )

# Keep protein ID and binary assignment
binary = binary[, c("protein", "above")]
myimmun_binary = binary

# Remove domains where the comparison could not be made
myimmun_binary = myimmun_binary %>%
  filter(!is.na(above))

# name
myimmun_binary$group = "Immune"


## regression
myimmun = domain_non_domain
imm_model  = lmer(Ratio ~ Ratio_non_domain + (1 | protein), data = domain_non_domain)

imm_int   = fixef(imm_model)[["(Intercept)"]]
imm_slope = fixef(imm_model)[["Ratio_non_domain"]]


immun_main = ggplot(
  domain_non_domain,
  aes(x = Ratio_non_domain, y = Ratio)
) +
  geom_abline(
    intercept = imm_int,
    slope = imm_slope,
    color = "#65c2a5",
    linewidth = 1
  ) +
  geom_point(size = 0.5, alpha = 1, color = "#65c2a5") +
  #geom_smooth(method = "lm", color = "#65c2a5") +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  theme_classic(base_size = 10) +
  labs(
    x = "R/C ratio (non-domain)",
    y = ""
  )+ coord_cartesian(xlim = c(0, 10), ylim = c(0, 10))#+ annotate( "text", x = 2.5, y = 7.5, label = paste0("%", round(above_line, 2)) )


rm(list = setdiff(ls(), c("Pleio_nondom","Pleio_nondom_inset" ,"Pleio_main","pabove","ptotal",
                          "immun_nondom","immun_nondom_inset" ,"immun_main","iabove","itotal",
                          "mypleio","myimmun",
                          "mypleio_binary","myimmun_binary")))










######################### 
#devel
#########################


setwd("..")
setwd("devel")


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


df = df %>%
  filter(protein != "non_domain_all")

non_domains = df %>%
  filter(domain %in% c("non_domain", "non_domain_all")) %>%
  select(protein, domain, Ratio) %>%
  pivot_wider(
    names_from = domain,
    values_from = Ratio
  ) %>%
  filter(
    is.finite(non_domain),
    is.finite(non_domain_all)
  )


# calculate deviation from y=x
non_domains= non_domains %>%
  mutate(
    deviation = non_domain - non_domain_all
  )

# Deviation histogram
deviation_hist = ggplot(non_domains, aes(x = deviation)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    bins = 100
  ) +
  # geom_vline(
  #   xintercept = 0,
  #   linetype = "dashed",
  #   linewidth = 0.3
  # ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1)
  ) +
  theme_classic(base_size = 7) +
  labs(
    x = "Deviation from Y = X",
    y = "Proportion"
  )


devel_nondom = ggplot(non_domains, aes(x = non_domain_all, y = non_domain)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  coord_equal() +
  theme_classic(base_size = 13) +
  labs(
    x = "R/C ratio (non-domain-all)",
    y = ""
  )+ggtitle('Developmental')

devel_nondom_inset = devel_nondom +
  inset_element(
    deviation_hist,
    left = 0.58,
    bottom = 0.55,
    right = 0.98,
    top = 0.98
  )


#### remove those with dS<3
mydnds = read.csv("/Users/danial/Documents/postdoc_research/Domains/My_approach/seq_3_groups/main/updated/AllDomains/sixspeices/a_a_changes/devel/devel_dnds.csv",row.names = 1)
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

#write
#write.csv(domain_non_domain, "/Users/danial/Documents/postdoc_research/Domains/My_approach/seq_3_groups/main/updated/AllDomains/sixspeices/a_a_changes/devel/devel_RC.csv")



# above y=x
dabove = nrow(domain_non_domain[domain_non_domain$Ratio > domain_non_domain$Ratio_non_domain, ])
dtotal = nrow(domain_non_domain)
above_line = dabove/dtotal* 100


# bindary assignment for each domain
binary = domain_non_domain %>%
  mutate(
    above = as.integer(Ratio > Ratio_non_domain)
  )

# Keep protein ID and binary assignment
binary = binary[, c("protein", "above")]
mydevel_binary = binary

# Remove domains where the comparison could not be made
mydevel_binary = mydevel_binary %>%
  filter(!is.na(above))

# name
mydevel_binary$group = "Developmental"


## regression
mydevel = domain_non_domain
dev_model  = lmer(Ratio ~ Ratio_non_domain + (1 | protein), data = domain_non_domain)

dev_int   = fixef(dev_model)[["(Intercept)"]]
dev_slope = fixef(dev_model)[["Ratio_non_domain"]]


devel_main = ggplot(
  domain_non_domain,
  aes(x = Ratio_non_domain, y = Ratio)
) +
  geom_abline(
    intercept = dev_int,
    slope = dev_slope,
    color = "#8da0cb",
    linewidth = 1
  ) +
  geom_point(size = 0.5, alpha = 1, color = "#8da0cb") +
  #geom_smooth(method = "lm", color = "#8da0cb") +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  theme_classic(base_size = 10) +
  labs(
    x = "R/C ratio (non-domain)",
    y = ""
  )+ coord_cartesian(xlim = c(0, 10), ylim = c(0, 10))#+ annotate( "text", x = 2.5, y = 7.5, label = paste0("%", round(above_line, 2)) )

rm(list = setdiff(ls(), c("Pleio_nondom", "Pleio_nondom_inset" ,"Pleio_main","pabove","ptotal",
                          "immun_nondom", "immun_nondom_inset" ,"immun_main","iabove","itotal",
                          "devel_nondom", "devel_nondom_inset" ,"devel_main","dabove","dtotal",
                          "mypleio","myimmun","mydevel",
                          "mypleio_binary","myimmun_binary","mydevel_binary")))









#################### 
#rando
####################


setwd('..')
setwd('rando')


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




non_domains = df %>%
  filter(domain %in% c("non_domain", "non_domain_all")) %>%
  select(protein, domain, Ratio) %>%
  pivot_wider(
    names_from = domain,
    values_from = Ratio
  )


# calculate deviation from y=x
non_domains= non_domains %>%
  mutate(
    deviation = non_domain - non_domain_all
  )

# Deviation histogram
deviation_hist = ggplot(non_domains, aes(x = deviation)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    bins = 100
  ) +
  # geom_vline(
  #   xintercept = 0,
  #   linetype = "dashed",
  #   linewidth = 0.3
  # ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1)
  ) +
  theme_classic(base_size = 7) +
  labs(
    x = "Deviation from Y = X",
    y = "Proportion"
  )


rando_nondom = ggplot(non_domains, aes(x = non_domain_all, y = non_domain)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  coord_equal() +
  theme_classic(base_size = 13) +
  labs(
    x = "R/C ratio (non-domain-all)",
    y = ""
  )+ggtitle('Random')

rando_nondom_inset = rando_nondom +
  inset_element(
    deviation_hist,
    left = 0.58,
    bottom = 0.55,
    right = 0.98,
    top = 0.98
  )

#### remove those with dS<3
mydnds = read.csv("/Users/danial/Documents/postdoc_research/Domains/My_approach/seq_3_groups/main/updated/AllDomains/sixspeices/a_a_changes/rando/rando_dnds.csv",row.names = 1)
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


#write
#write.csv(domain_non_domain, "/Users/danial/Documents/postdoc_research/Domains/My_approach/seq_3_groups/main/updated/AllDomains/sixspeices/a_a_changes/rando/rando_RC.csv")


# above y=x
rabove = nrow(domain_non_domain[domain_non_domain$Ratio > domain_non_domain$Ratio_non_domain, ])
rtotal = nrow(domain_non_domain)
above_line = rabove/rtotal* 100

# bindary assignment for each domain
binary = domain_non_domain %>%
  mutate(
    above = as.integer(Ratio > Ratio_non_domain)
  )

# Keep protein ID and binary assignment
binary = binary[, c("protein", "above")]
myrando_binary = binary

# Remove domains where the comparison could not be made
myrando_binary = myrando_binary %>%
  filter(!is.na(above))

# name
myrando_binary$group = "Random"

## regression
myrando = domain_non_domain
ran_model  = lmer(Ratio ~ Ratio_non_domain + (1 | protein), data = domain_non_domain)

ran_int   = fixef(ran_model)[["(Intercept)"]]
ran_slope = fixef(ran_model)[["Ratio_non_domain"]]


rando_main = ggplot(
  domain_non_domain,
  aes(x = Ratio_non_domain, y = Ratio)
) +
  geom_abline(
    intercept = ran_int,
    slope = ran_slope,
    color = "grey",
    linewidth = 1
  ) +
  geom_point(size = 0.5, alpha = 1, color = "grey") +
  #geom_smooth(method = "lm", color = "grey") +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  theme_classic(base_size = 10) +
  labs(
    x = "R/C ratio (non-domain)",
    y = ""
  )+ coord_cartesian(xlim = c(0, 10), ylim = c(0, 10))#+ annotate( "text", x = 2.5, y = 7.5, label = paste0("%", round(above_line, 2)) )



rm(list = setdiff(ls(), c("Pleio_nondom", "Pleio_nondom_inset","Pleio_main","pabove","ptotal",
                          "immun_nondom", "immun_nondom_inset","immun_main","iabove","itotal",
                          "devel_nondom", "devel_nondom_inset","devel_main","dabove","dtotal",
                          "rando_nondom", "rando_nondom_inset","rando_main","rabove","rtotal",
                          "mypleio","myimmun","mydevel","myrando",
                          "mypleio_binary","myimmun_binary","mydevel_binary","myrando_binary")))
# 4x16 pdf
(   Pleio_nondom_inset |
    immun_nondom_inset |
    devel_nondom_inset |
    rando_nondom_inset)

grid.arrange(Pleio_main,
             immun_main,
             devel_main,
             rando_main, nrow =1)







######################################## regressions

length(unique(mydevel$protein))
length(unique(myimmun$protein))
length(unique(mypleio$protein))
length(unique(myrando$protein))


dev_model  = lmer(Ratio ~ Ratio_non_domain + (1 | protein), data = mydevel)
imm_model  = lmer(Ratio ~ Ratio_non_domain + (1 | protein), data = myimmun)
ple_model  = lmer(Ratio ~ Ratio_non_domain + (1 | protein), data = mypleio)
ran_model  = lmer(Ratio ~ Ratio_non_domain + (1 | protein), data = myrando)

mydevel$group = "Developmental"
myimmun$group = "Immune"
mypleio$group = "Pleiotropic"
myrando$group = "Random"

all_dat = rbind(
  mydevel,
  myimmun,
  mypleio,
  myrando
)


#no evidence of different slopes among the groups.
combined_model = lmerTest::lmer(
  Ratio ~ Ratio_non_domain * group + (1 | protein),
  data = all_dat
)

summary(combined_model)
anova(combined_model)
emtrends(combined_model, pairwise ~ group, var = "Ratio_non_domain")


########################## separate linear mixed model fits
##########################################################

ple_int   = fixef(ple_model)[["(Intercept)"]]
ple_slope = fixef(ple_model)[["Ratio_non_domain"]]
imm_int   = fixef(imm_model)[["(Intercept)"]]
imm_slope = fixef(imm_model)[["Ratio_non_domain"]]
dev_int   = fixef(dev_model)[["(Intercept)"]]
dev_slope = fixef(dev_model)[["Ratio_non_domain"]]
ran_int   = fixef(ran_model)[["(Intercept)"]]
ran_slope = fixef(ran_model)[["Ratio_non_domain"]]

round(range(c(ple_slope, imm_slope, dev_slope, ran_slope)),2)

lines_lmxmodel <- data.frame(
  Group = c("Developmental", "Immunity", "Pleiotropy", "Random"),
  intercept = c(
    dev_int,
    imm_int,
    ple_int,
    ran_int
    
  ),
  slope = c(
    dev_slope,
    imm_slope,
    ple_slope,
    ran_slope
    
  )
)

lines_lmxmodel$Group = factor(lines_lmxmodel$Group, c('Random', 'Immunity', 'Developmental', 'Pleiotropy') )

Regline_plot = ggplot(lines_lmxmodel) +
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
    x = "R/C (non-domain)",
    y = "Predicted R/C (domain)",
    color = NULL
  ) +
  theme(
    legend.position = "",
  )+
  scale_x_continuous(
    breaks = seq(0, 10.5, 2.5),
    limits = c(0, 10.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 10.5, 2.5),
    limits = c(0, 10.5),
    expand = c(0, 0)
  )+geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) 

grid.arrange(Pleio_main,
             immun_main,
             devel_main,
             rando_main, 
             Regline_plot,nrow =1)
