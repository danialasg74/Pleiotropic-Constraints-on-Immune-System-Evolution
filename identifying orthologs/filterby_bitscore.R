setwd('/Users/danial/Documents/postdoc_research/Domains/My_approach/seq_3_groups/Finding_orthologs/D_simulans')

pleio = read.delim("pleio.txt", header = FALSE)
pleio_order = pleio[order(-pleio$V12), ]
pleio_order_uniq = subset(pleio_order, !duplicated(pleio_order$V1))
pleio_quaility = subset(pleio_order_uniq, pleio_order_uniq$V12>100)
sum(pleio_quaility$V11>0.05)

immun = read.delim("immun.txt", header = FALSE)
immun_order = immun[order(-immun$V12), ]
immun_order_uniq = subset(immun_order, !duplicated(immun_order$V1))
immun_quaility = subset(immun_order_uniq, immun_order_uniq$V12>100)
sum(immun_quaility$V11>0.05)

devel = read.delim("devel.txt", header = FALSE)
devel_order = devel[order(-devel$V12), ]
devel_order_uniq = subset(devel_order, !duplicated(devel_order$V1))
devel_quaility = subset(devel_order_uniq, devel_order_uniq$V12>100)
sum(devel_quaility$V11>0.05)


random = read.delim("random.txt", header = FALSE)
random_order = random[order(-random$V12), ]
random_order_uniq = subset(random_order, !duplicated(random_order$V1))
random_quaility = subset(random_order_uniq, random_order_uniq$V12>100)
sum(random_quaility$V11>0.05)





write.csv(pleio_quaility,'pleio_match.csv')
write.csv(immun_quaility,'immun_match.csv')
write.csv(devel_quaility,'devel_match.csv')
write.csv(random_quaility,'rando_match.csv')


