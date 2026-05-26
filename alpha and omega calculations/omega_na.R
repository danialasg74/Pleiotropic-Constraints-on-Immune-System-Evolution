# Base path (up to RAL)

############ go to the population directory

base_dir <- "../MUN"

############ each folder contains MultiDFE outputs

groups <- c("pleio", "immun", "devel", 'random')

# ---- Helper function ----
get_fix_prob <- function(f, group_name) {
  line <- readLines(f, warn = FALSE)
  fields <- strsplit(line, "\t")[[1]]
  kv <- strsplit(fields, ":")
  
  keys <- vapply(kv, `[`, character(1), 1)
  vals <- vapply(kv, `[`, character(1), 2)
  
  fp <- vals[keys == "fix_prob"]
  fp <- if (length(fp) == 0) NA_real_ else as.numeric(fp)
  
  boot <- suppressWarnings(
    as.integer(sub(".*bootstrap_([0-9]+).*", "\\1", basename(f)))
  )
  
  data.frame(
    group = group_name,
    file = basename(f),
    boot = boot,
    fix_prob = fp,
    stringsAsFactors = FALSE
  )
}

# ---- Loop over groups ----
all_fix_probs <- do.call(
  rbind,
  lapply(groups, function(g) {
    
    dir_path <- file.path(base_dir, g, "bootstrap")
    files <- list.files(dir_path, pattern = "\\.MAXL\\.out$", full.names = TRUE)
    
    do.call(
      rbind,
      lapply(files, get_fix_prob, group_name = g)
    )
  })
)

# Sort 
all_fix_probs <- all_fix_probs[order(all_fix_probs$group, all_fix_probs$boot), ]

pleio_med = median(subset(all_fix_probs, all_fix_probs$group == 'pleio')$fix_prob)
immun_med = median(subset(all_fix_probs, all_fix_probs$group == 'immun')$fix_prob)
devel_med = median(subset(all_fix_probs, all_fix_probs$group == 'devel')$fix_prob)
rando_med = median(subset(all_fix_probs, all_fix_probs$group == 'random')$fix_prob)


############ plot ############

library(ggplot2)

all_fix_probs$group = factor(all_fix_probs$group, levels = c("random",'immun','pleio','devel'))

ggplot(all_fix_probs, aes(x = group, y = fix_prob, fill = group)) +
  geom_violin(trim = TRUE, alpha = 1) +
  geom_boxplot(width = 0.1,outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  scale_fill_manual(
    values = c(
      immun = "#65c2a5",
      pleio = "#fc8d62",
      devel = "#8da0cb",
      random = "grey"
    ),
    labels = c(
      immun = "Immunity",
      pleio = "Pleiotropic",
      devel = "Developmental",
      random = "Random"
    )
  )  +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "omega_na")+ylim(0,1)+scale_x_discrete(labels = c(
    immun = "Immunity",
    pleio = "Pleiotropic",
    devel = "Developmental",
    random = "Random"
  ))+annotate(
    "segment",
    x =1.4, xend = 2.5,
    y = immun_med, yend = immun_med,
    colour = "black", # pick one of your fill colors
    linetype = "dashed",
    linewidth = 0.3
  )+annotate(
    "segment",
    x = 2.3, xend = 3.4,
    y = pleio_med, yend = pleio_med,
    colour = "black", # pick one of your fill colors
    linetype = "dashed",
    linewidth = 0.3
  )+annotate(
    "segment",
    x = 3.2, xend = 4.3,
    y = devel_med, yend = devel_med,
    colour = "black", # pick one of your fill colors
    linetype = "dashed",
    linewidth = 0.3)+
  annotate(
    "segment",
    x =0.5, xend = 1.6,
    y = rando_med, yend = rando_med,
    colour = "black", # pick one of your fill colors
    linetype = "dashed",
    linewidth = 0.3)+ylim(0,0.1)+ annotate(
      "text",
      x = 1,
      y = 0.1,
      label = "a",
      colour = "black",
      size = 5
    )+
  annotate(
    "text",
    x = 2,
    y = 0.1,
    label = "b",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 3,
    y = 0.1,
    label = "c",
    colour = "black",
    size = 5
  )+
  annotate(
    "text",
    x = 4,
    y = 0.1,
    label = "d",
    colour = "black",
    size = 5)+scale_y_continuous(
      breaks = c(0, 0.025, 0.05, 0.075, 0.1),
      limits = c(0, 0.13)
    )



#immun → A
#pleio → B
#devel → C
#random → D

####### stat test
x = c(all_fix_probs$fix_prob)
grp = c(all_fix_probs$group)
kruskal.test(x ~ grp)
pairwise.wilcox.test(x, grp, p.adjust.method = "BH", exact = FALSE, na.action = na.omit)




########################## save MultiDFE-outputted probabilities ##########################


pleio = subset(all_fix_probs, all_fix_probs$group == 'pleio')[c(2,4)]
pleio$file = sub("\\.sfs\\.MAXL\\.out$", "", pleio$file)

write.table(pleio,
            file = "pleio.txt",
            sep = "\t",        # tab-separated (change to "," if you want CSV)
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)




immun = subset(all_fix_probs, all_fix_probs$group == 'immun')[c(2,4)]
immun$file = sub("\\.sfs\\.MAXL\\.out$", "", immun$file)

write.table(immun,
            file = "immun.txt",
            sep = "\t",        # tab-separated (change to "," if you want CSV)
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)




devel = subset(all_fix_probs, all_fix_probs$group == 'devel')[c(2,4)]
devel$file = sub("\\.sfs\\.MAXL\\.out$", "", devel$file)

write.table(devel,
            file = "devel.txt",
            sep = "\t",        # tab-separated (change to "," if you want CSV)
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)


random = subset(all_fix_probs, all_fix_probs$group == 'random')[c(2,4)]
random$file = sub("\\.sfs\\.MAXL\\.out$", "", random$file)

write.table(random,
            file = "random.txt",
            sep = "\t",        # tab-separated (change to "," if you want CSV)
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)