########################## Import data from Williams et al. 2023 supplemental data ##########################

library(readxl)
setwd('/Users/danial/Documents/postdoc_research/Domains/My_approach/seq_3_groups/Finding_ncbi_seq')
FlyIDs = read_excel(
  "17Jan2023_Supplemental_File_1.xlsx",
  sheet = "Supplemental_Table_1"
)


Pleiotropic = subset(FlyIDs, FlyIDs$Plei_immuneResponse_allDev=='Pleiotropic')
Immunitiees = subset(FlyIDs, FlyIDs$Plei_immuneResponse_allDev=='Immune_Non_Pleiotropic')
Development = subset(FlyIDs, FlyIDs$Plei_immuneResponse_allDev=='Developmental_Non_Pleiotropic')


#write.csv(Pleiotropic, 'Pleiotropic.csv')
#write.csv(Immunitiees, 'Imunity.csv')
#write.csv(Development, 'Development.csv')


########################## Find NCBI ID ##########################

rm(FlyIDs)
Pleiotropic = Pleiotropic[,1]
Immunitiees = Immunitiees[,1]
Development = Development[,1]

# import melanogaster annotation from biomart (https://metazoa.ensembl.org/biomart/martview/43bc13deb09cd272ab2a9a8db873332e)
mart = read.delim('/Users/danial/Documents/postdoc_research/Domains/My_approach/seq_3_groups/Finding_ncbi_seq/mart_export.txt')

#> head(mart)
#Gene.stable.ID Protein.stable.ID NCBI.gene..formerly.Entrezgene..ID NCBI.gene..formerly.Entrezgene..accession
#1    FBti0019141                                                   NA                                        NA
#2    FBti0019141                                                   NA                                        NA
#3    FBti0019141                                                   NA                                        NA
#4    FBgn0053832       FBpp0091078                            3772632                                   3772632
#5    FBgn0053832       FBpp0091078                            3772345                                   3772345
#6    FBgn0053832       FBpp0091078                            3771783                                   3771783


names(mart)[1] = names(Pleiotropic)[1]

Pleiotropic_mart = merge(Pleiotropic, mart)
Immunitiees_mart = merge(Immunitiees, mart)
Development_mart = merge(Development, mart)


Pleiotropic_mart = subset(Pleiotropic_mart, !duplicated(Pleiotropic_mart$FlyBaseID))
Immunitiees_mart = subset(Immunitiees_mart, !duplicated(Immunitiees_mart$FlyBaseID))
Development_mart = subset(Development_mart, !duplicated(Development_mart$FlyBaseID))



########################## which genes have no correspondent Id ##########################


setdiff(Pleiotropic$FlyBaseID, Pleiotropic_mart$FlyBaseID)
setdiff(Immunitiees$FlyBaseID, Immunitiees_mart$FlyBaseID)
setdiff(Development$FlyBaseID, Development_mart$FlyBaseID)

paste0('Pleiotropic:', length(setdiff(Pleiotropic$FlyBaseID, Pleiotropic_mart$FlyBaseID)))
paste0('Immunitiees:', length(setdiff(Immunitiees$FlyBaseID, Immunitiees_mart$FlyBaseID)))
paste0('Development:', length(setdiff(Development$FlyBaseID, Development_mart$FlyBaseID)))


########################## NCBI annotation ##########################

Pleiotropic_mart = Pleiotropic_mart[,-4]
Immunitiees_mart = Immunitiees_mart[,-4]
Development_mart = Development_mart[,-4]

NCBI = read.delim('ncbi_dataset.tsv')
NCBI = NCBI[,c(7,8,12,14,15)]

names(Pleiotropic_mart)[3] = names(NCBI)[2]
names(Immunitiees_mart)[3] = names(NCBI)[2]
names(Development_mart)[3] = names(NCBI)[2]

Pleiotropic_ncbi = merge(Pleiotropic_mart, NCBI, by = 'Gene.ID')
Immunitiees_ncbi = merge(Immunitiees_mart, NCBI, by = 'Gene.ID')
Development_ncbi = merge(Development_mart, NCBI, by = 'Gene.ID')


Pleiotropic_ncbi_long = Pleiotropic_ncbi[order(-Pleiotropic_ncbi$Protein.length), ]
Immunitiees_ncbi_long = Immunitiees_ncbi[order(-Immunitiees_ncbi$Protein.length), ]
Development_ncbi_long = Development_ncbi[order(-Development_ncbi$Protein.length), ]


Pleiotropic_ncbi_long = subset(Pleiotropic_ncbi_long, !duplicated(Pleiotropic_ncbi_long$Gene.ID))
Immunitiees_ncbi_long = subset(Immunitiees_ncbi_long, !duplicated(Immunitiees_ncbi_long$Gene.ID))
Development_ncbi_long = subset(Development_ncbi_long, !duplicated(Development_ncbi_long$Gene.ID))


write.csv(Pleiotropic_ncbi_long, 'Pleiotropic_ncbi.csv')
write.csv(Immunitiees_ncbi_long, 'Immunity_ncbi.csv')
write.csv(Development_ncbi_long, 'Development_ncbi.csv')
