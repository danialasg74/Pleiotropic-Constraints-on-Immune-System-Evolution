library(jsonlite)
protein_folders <- list.dirs(".", recursive = FALSE)

# Function to compute N / (len/3) for ANY element (protein or domain)
compute_ratio <- function(name, change_list, length_list) {
  # bail out safely if a file was missing/filtered out
  if (is.null(change_list[[name]]) || is.null(length_list[[name]])) {
    return(NA_real_)
  }
  df <- data.frame(change_list[[name]])
  N  <- sum(df[, ncol(df)] < 0)
  
  len <- as.numeric(data.frame(length_list[[name]])$V1)
  len_nt <- len / 3  # convert nt → aa length
  
  N / len_nt
}

for (pf in protein_folders) {
  message("Processing: ", pf)
  
  changes <- list.files(
    path = pf,
    pattern = "^branchchanges_all\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  LENGTH <- list.files(
    path = pf,
    pattern = "length",
    full.names = TRUE,
    recursive = TRUE
  )
  
  # --- NEW: skip empty files to avoid read.* errors ---
  if (length(changes)) {
    changes <- changes[file.info(changes)$size > 0]
  }
  if (length(LENGTH)) {
    LENGTH  <- LENGTH[file.info(LENGTH)$size > 0]
  }
  # ----------------------------------------------------
  
  # Read all CSVs found under this protein folder
  change_list <- lapply(changes, read.csv)
  length_list <- lapply(LENGTH,  read.table)
  
  # Name each list element
  names(change_list) <- vapply(changes, function(f) {
    f_norm  <- gsub("\\\\", "/", f)
    pf_norm <- gsub("\\\\", "/", pf)
    parent_dir <- gsub("\\\\", "/", dirname(f))
    if (normalizePath(parent_dir, winslash = "/") ==
        normalizePath(pf_norm,  winslash = "/")) {
      basename(pf_norm)
    } else {
      basename(parent_dir)
    }
  }, character(1))
  
  names(length_list) <- vapply(LENGTH, function(f) {
    f_norm  <- gsub("\\\\", "/", f)
    pf_norm <- gsub("\\\\", "/", pf)
    parent_dir <- gsub("\\\\", "/", dirname(f))
    if (normalizePath(parent_dir, winslash = "/") ==
        normalizePath(pf_norm,  winslash = "/")) {
      basename(pf_norm)
    } else {
      basename(parent_dir)
    }
  }, character(1))
  
  protein_name <- basename(pf)
  domain_names <- setdiff(names(change_list), protein_name)
  
  # --- Protein ---
  prot_ratio <- compute_ratio(protein_name, change_list, length_list)
  
  # --- Domains ---
  domain_ratios <- sapply(domain_names, compute_ratio,
                          change_list = change_list,
                          length_list = length_list)
  
  result <- data.frame(
    name  = c(protein_name, domain_names),
    ratio = c(prot_ratio, as.numeric(domain_ratios))
  )
  write_json(result, file.path(pf, "normalized_radical_all.json"), row.names = FALSE)
}
