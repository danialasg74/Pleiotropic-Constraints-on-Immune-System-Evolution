#load packages
library(dplyr)
#get all subfolder names
dirs = list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
for(dir_indx in 1:length(dirs)){
  #choose a subfolder
  currnt_dir = dirs[dir_indx]
  #import domain
  if (file.exists(paste0(currnt_dir, '/domains.tsv'))) {
    domain = read.delim(paste0(currnt_dir, '/domains.tsv'))
  }
  
  #remove short (i.e., < 60 a.a domains)
  domain = domain[abs(domain$ali_from - domain$ali_to) > 20 ,]
  
  #remove uncertain domains (i.e., E val > 10^-5)
  domain = domain[domain$i_evalue <  10^-5,]
  
  
  domain <- domain %>%
    group_by(query, pfam_acc) %>%
    mutate(pfam_acc = if (n() > 1) paste0(pfam_acc, "_", row_number()) else pfam_acc) %>%
    ungroup()
  
  #write the unique domain
  write.table(domain, paste0(currnt_dir,"/domains_20.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

}




