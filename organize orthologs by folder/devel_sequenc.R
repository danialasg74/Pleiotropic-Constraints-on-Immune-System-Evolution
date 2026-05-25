############### Dmel Ids and seqs for devele genes ###############
library(Biostrings)

dmel_devel_Id = read.csv('Development_ncbi.csv')
dmel_devel_Id =  dmel_devel_Id$Protein.accession[dmel_devel_Id$Protein.accession != ""]

mel_cds = readDNAStringSet("dmel_devel.fasta")

###############  Ids 5 other species; developmental genes ###############

setwd('../D_ananassae')
ana_devel_Id = read.csv('devel_match.csv')
ana_cds = readDNAStringSet("ana_devel.fasta")

setwd('../D_erecta')
ere_devel_Id = read.csv('devel_match.csv')
ere_cds = readDNAStringSet("ere_devel.fasta")

setwd('../D_sechellia')
sec_devel_Id = read.csv('devel_match.csv')
sec_cds = readDNAStringSet("sec_devel.fasta")

setwd('../D_simulans')
sim_devel_Id = read.csv('devel_match.csv')
sim_cds = readDNAStringSet("sim_devel.fasta")

setwd('../D_yakuba')
yak_devel_Id = read.csv('devel_match.csv')
yak_cds = readDNAStringSet("yak_devel.fasta")

###############  find orthologs ###############


##### functions 

### get the ortholog ###
getV2 = function(df, key) {
  hits <- df$V2[df$V1 == key]
  if (length(hits) == 0) NA_character_ else hits[1]  # take first if duplicates
}

#Get the sequence for a unique Id
getsequence = function(cds, id) {
  if (length(id) == 0 || is.na(id) || id == "") return(cds[0])
  cds[grepl(id, names(cds), fixed = TRUE)]
}


for(i in 1:length(dmel_devel_Id)){
  

  ######### i = the ith ID from deml #########
  #choose a Dmel gene
  key = dmel_devel_Id[i]
  
  #find orthologs
  IDs = data.frame(
    mel = key,
    ana = getV2(ana_devel_Id, key),
    ere = getV2(ere_devel_Id, key),
    sec = getV2(sec_devel_Id, key),
    sim = getV2(sim_devel_Id, key),
    yak = getV2(yak_devel_Id, key),
    stringsAsFactors = FALSE
  )
  
  ###############  find sequences ###############
  
  
  
  fasta_mel = getsequence(mel_cds, IDs$mel)
  fasta_ana = getsequence(ana_cds, IDs$ana)
  fasta_ere = getsequence(ere_cds, IDs$ere)
  fasta_sec = getsequence(sec_cds, IDs$sec)
  fasta_sim = getsequence(sim_cds, IDs$sim)
  fasta_yak = getsequence(yak_cds, IDs$yak)
  
  if (length(names(fasta_mel)) > 0) names(fasta_mel) = "dmel"
  if (length(names(fasta_ana)) > 0) names(fasta_ana) = "dana"
  if (length(names(fasta_ere)) > 0) names(fasta_ere) = "dere"
  if (length(names(fasta_sec)) > 0) names(fasta_sec) = "dsec"
  if (length(names(fasta_sim)) > 0) names(fasta_sim) = "dsim"
  if (length(names(fasta_yak)) > 0) names(fasta_yak) = "dyak"
  
  
  fasta_list = list(fasta_mel,
                    fasta_ana,
                    fasta_ere,
                    fasta_sec,
                    fasta_sim,
                    fasta_yak)
  
  combined = do.call(c, fasta_list)
  
  dir.create(paste(key))
  setwd(paste(key))
  writeXStringSet(combined, "seq.fasta")
}












