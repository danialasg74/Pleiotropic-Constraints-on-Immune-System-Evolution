library(Biostrings)

#Protein folders
dirs <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)

base_dir <- getwd()  # remember where we started

for (dir_indx in seq_along(dirs)) {
  cat("Processing:", dir_indx, dirs[dir_indx], "\n")
  
  # go to subfolder (i.e., protein)
  currnt_dir <- dirs[dir_indx]
  setwd(file.path(base_dir, currnt_dir))
  
  # init placeholders to avoid "object not found"
  protein <- NULL
  domain_sub <- NULL
  
  # get the protein sequences for all species and nodes
  if (file.exists('pro.fasta')) {
    protein <- readAAStringSet('pro.fasta')
  }
  
  # get domain coordinates
  if (file.exists('domains_20.tsv')) {
    domain_sub <- read.delim('domains_20.tsv')
  }
  
  if (!is.null(protein) && length(protein) >= 1 &&
      !is.null(domain_sub) && nrow(domain_sub) >= 1) {
    
    
    # for every domain
    fasta_lines <- character(0)
    
    for (i in 1:length(unique(domain_sub$query))) {
      # name of the domain
      myquery = unique(domain_sub$query)[i]
      
      if (!myquery %in% names(protein)) next
      
      # get domain coords for all species
      domain_species = subset(domain_sub, query == myquery)
    
      
      # Make domain intervals
      domain_ranges = IRanges(
        start = domain_species$ali_from,
        end   = domain_species$ali_to
      )
      
      # Merge overlapping domains
      domain_ranges =  reduce(domain_ranges)
      
      # Get regions outside the domains
      non_domain_ranges = gaps(
        domain_ranges,
        start = 1,
        end = width(protein[myquery])
      )
      
      # Extract the non-domain sequences
      non_domain_views = Views(
        protein[myquery][[1]],
        non_domain_ranges
      )
      
      
      if (length(non_domain_views) == 0) next
      # Concatenate all non-domain regions
      my_concatenate = do.call(xscat, as.list(non_domain_views))
      
      
      
      # Create FASTA header
      header = paste0(">", myquery)
      
      # Add header and sequence to FASTA output
      fasta_lines = c(fasta_lines, header, as.character(my_concatenate))
      
    }
    
    # only write if we actually collected something
    if (length(fasta_lines) > 0) {
      #create a non-domain directory
      dir.create('non_domain', showWarnings = FALSE)
      setwd('non_domain')
      #write the fasta
      writeLines(fasta_lines, "non_domain_pro.fasta")
    }
    
    # go back to the base directory
    setwd(base_dir)
  }
}





