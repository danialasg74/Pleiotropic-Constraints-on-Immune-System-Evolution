library(Biostrings)
#domain_name
get_subseq <- function(query_name, domain_tbl, codon_list) {
  ali_from <- domain_tbl$ali_from[domain_tbl$query == query_name]
  ali_to   <- domain_tbl$ali_to[domain_tbl$query == query_name]
  subseq(codon_list[[query_name]], start = ali_from, end = ali_to)
}

dirs <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)

base_dir <- getwd()  # remember where we started

for (dir_indx in seq_along(dirs)) {
  cat("Processing:", dir_indx, dirs[dir_indx], "\n")
  
  # go to subfolder (i.e., protein)
  currnt_dir <- dirs[dir_indx]
  setwd(file.path(base_dir, currnt_dir))
  
  # init placeholders to avoid "object not found"
  codons <- NULL
  domain_sub <- NULL
  
  # get the protein sequences for all species and nodes
  if (file.exists('seq.fasta')) {
    codons <- readDNAStringSet('seq.fasta')
  }
  
  # get domain coordinates
  if (file.exists('domains_20.tsv')) {
    domain_sub <- read.delim('domains_20.tsv')
  }
  
  if (!is.null(codons) && length(codons) >= 1 &&
      !is.null(domain_sub) && nrow(domain_sub) >= 1) {
    
    pfams <- unique(domain_sub$pfam_acc)
    
    # for every domain
    for (i in seq_along(pfams)) {
      # name of the domain
      domain_name <- pfams[i]
      
      # get domain coords for all species
      domain_tbl <- subset(domain_sub, pfam_acc == domain_name)
      
      #converting aa coordinates to nucleotides
      domain_tbl$ali_from = (domain_tbl$ali_from - 1)*3 + 1
      domain_tbl$ali_to   = domain_tbl$ali_to*3
      
      
      # skip if empty for safety
      if (nrow(domain_tbl) == 0) next
      
      # make a subfolder using the domain name
      dir.create(paste0(domain_name,'_20'), showWarnings = FALSE)
      
      # go to the domain folder
      setwd(paste0(domain_name,'_20'))
      
      # one FASTA per domain
      fasta_lines <- character(0)
      
      for (j in seq_len(nrow(domain_tbl))) {
        qname <- domain_tbl$query[j]
        
        # ensure the sequence exists in the AAStringSet
        if (!qname %in% names(codons)) next
        
        seq_j  <- get_subseq(qname, domain_tbl, codons)
        header <- paste0(">", qname)
        fasta_lines <- c(fasta_lines, header, as.character(seq_j))
      }
      
      # only write if we actually collected something
      if (length(fasta_lines) > 0) {
        writeLines(fasta_lines, paste0(domain_name, ".fasta"))
      }
      
      # go back to the protein folder
      setwd("..")
    }
  }
  
  # always go back to the base directory before next iteration
  setwd(base_dir)
}
