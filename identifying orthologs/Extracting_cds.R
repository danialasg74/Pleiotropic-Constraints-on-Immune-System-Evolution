setwd('../D_simulans')
library(Biostrings)
fasta = readDNAStringSet("cds_Dsim.fna")
fasta[[1]]
names(fasta)[1]

pleio = read.csv('pleio_match.csv',row.names = 1) 
immun = read.csv('immun_match.csv',row.names = 1)
devel = read.csv('devel_match.csv',row.names = 1)
rando = read.csv('rando_match.csv',row.names = 1)

pleio_seq = fasta[grepl(paste(pleio$V2, collapse = "|"), names(fasta))]
immun_seq = fasta[grepl(paste(immun$V2, collapse = "|"), names(fasta))]
devel_seq = fasta[grepl(paste(devel$V2, collapse = "|"), names(fasta))]
rando_seq = fasta[grepl(paste(rando$V2, collapse = "|"), names(fasta))]




writeXStringSet(pleio_seq, filepath = "sim_pleio.fasta")
writeXStringSet(immun_seq, filepath = "sim_immun.fasta")
writeXStringSet(devel_seq, filepath = "sim_devel.fasta")
writeXStringSet(rando_seq, filepath = "sim_rando.fasta")
