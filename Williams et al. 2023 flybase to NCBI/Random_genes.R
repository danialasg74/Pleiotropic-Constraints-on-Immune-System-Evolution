################## 3 classes of genes ##################

setwd('...')
devel = read.csv('Non_Pleiotropic_Dev.csv')
immun = read.csv('Non_Pleiotropic_immune.csv')
pleio = read.csv('Pleiotropic.csv')


##total number of genes
SUM_three_rows  = nrow(devel) + nrow(immun) + nrow(pleio)
##average number of genes
Average_numb = SUM_three_rows/3
Average_numb = round(Average_numb)

################## all genes ##################


##import all genes and annotation
All_genes = read.delim('../fbgn_ncbi.txt')
All_gene_annot = read.delim('../ncbi_dataset.tsv')

## remove MT chromosomes
All_gene_annot = All_gene_annot[All_gene_annot$Chromosome != "MT", ]


##annotate genes
names(All_genes)[2] = names(All_gene_annot)[8]
All_genes_final = merge(All_genes, All_gene_annot, by =names(All_gene_annot)[8]  )

#choose protein coding and unique
All_genes_prot = subset(All_genes_final, All_genes_final$Gene.Type == "protein-coding")

All_genes_prot_unique = subset(All_genes_prot,!duplicated(All_genes_prot$Gene.stable.ID))

################## 3 classes combined ##################

Threecombined = c(devel$FlyBaseID,
  immun$FlyBaseID,
  pleio$FlyBaseID)


################## 3 classes excluded from "All_genes" ##################

not_in_three = setdiff(All_genes_prot_unique$Gene.stable.ID, Threecombined)
not_in_three = unique(not_in_three)

## get the N (average of three classes) number of genes randomly from not_in_three
set.seed(1)  # for reproducibility
random_genes = sample(not_in_three, Average_numb)

################## find prot id for random genes ##################

## annotate
random_genes_annot = All_genes_final[All_genes_final$Gene.stable.ID %in% random_genes,]

## sort by prot length
random_genes_annot_sorted = random_genes_annot[order(-random_genes_annot$Protein.length), ]

## get unique fbgn ids which yields the largest prot id too
random_genes_annot_unique = subset(random_genes_annot, !duplicated(random_genes_annot$Gene.stable.ID))

##proportion of chromosomes in random_genes_annot_unique

chromosome_vec = c()
number_vec = c()

for(i in 1:length(unique(random_genes_annot_unique$Chromosome))){
chromosome = unique(random_genes_annot_unique$Chromosome)[i]
number = nrow(subset(random_genes_annot_unique, random_genes_annot_unique$Chromosome == unique(random_genes_annot_unique$Chromosome)[i]))

chromosome_vec = c(chromosome_vec, chromosome)
number_vec = c(number_vec, number)
}

Num_Chro = data.frame(chrom = chromosome_vec, number = number_vec)
Num_Chro

barplot(
  Num_Chro$number,
  names.arg = Num_Chro$chrom,
  main = "#Random Genes",
  xlab = "Chromosome",
  ylab = "Number of Genes",
  ylim = c(0,300)
)

### in the ncbi annotation we have: 

## > nrow(All_gene_annot[All_gene_annot$Chromosome == "X", ])
## [1] 6045
## > nrow(All_gene_annot[All_gene_annot$Chromosome == "2L", ])
## [1] 6667
## > nrow(All_gene_annot[All_gene_annot$Chromosome == "2R", ])
## [1] 7010
## > nrow(All_gene_annot[All_gene_annot$Chromosome == "3L", ])
## [1] 6763
## > nrow(All_gene_annot[All_gene_annot$Chromosome == "3R", ])
## [1] 8175
## > nrow(All_gene_annot[All_gene_annot$Chromosome == "4", ])
## [1] 341
## > nrow(All_gene_annot[All_gene_annot$Chromosome == "Y", ])
## [1] 125



################## Export ##################

write.csv(random_genes_annot_unique, "random_genes.csv")


