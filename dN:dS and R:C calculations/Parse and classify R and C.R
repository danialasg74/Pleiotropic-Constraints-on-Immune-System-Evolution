# ---- Parse branch-change file and extract coordinates by species ----
library(jsonlite)
parse_branch_changes <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
  # indices of "Branch N:" lines
  b_idx <- grep("^Branch\\s+\\d+:", lines)
  if (length(b_idx) == 0) stop("No 'Branch' headers found.")
  
  # branch IDs
  b_id <- as.integer(sub("^Branch\\s+(\\d+):.*$", "\\1", lines[b_idx]))
  
  # species code in header if present (e.g., '(dmel)' right after the range)
  # We'll capture the first parenthetical group after the range "##..##" and before (n= ...),
  # which is where species tags appear in your example.
  get_species <- function(s) {
    # non-greedy up to first (...) and pull what's inside
    m <- regexpr("^Branch\\s+\\d+:.*?\\(([^)]+)\\)", s, perl = TRUE)
    if (m[1] == -1) return(NA_character_)
    val <- regmatches(s, m)
    sp  <- sub("^.*?\\(([^)]+)\\).*?$", "\\1", val, perl = TRUE)
    # If this first (...) is actually "n=..., s=..." then discard
    if (grepl("^\\s*n\\s*=|^\\s*s\\s*=", sp)) return(NA_character_)
    sp
  }
  species <- vapply(lines[b_idx], get_species, character(1))
  
  # end index for each branch block (start of next branch or end of file + 1)
  next_idx <- c(b_idx[-1], length(lines) + 1)
  
  # helper: extract change lines and fields from a block
  # Example change line:
  # "  400 GAC (D 0.999) -> GAG (E)"
  change_re <- "^\\s*(\\d+)\\s+([ACGT]{3})\\s+\\(([^)]+)\\)\\s*->\\s*([ACGT]{3})\\s+\\(([^)]+)\\)"
  
  out <- lapply(seq_along(b_idx), function(i) {
    block <- lines[(b_idx[i] + 1):(next_idx[i] - 1)]
    if (!length(block)) return(NULL)
    # keep only lines that look like changes
    keep <- grep(change_re, block)
    if (!length(keep)) return(NULL)
    cl <- block[keep]
    
    coord <- as.integer(sub(change_re, "\\1", cl, perl = TRUE))
    from_codon <- sub(change_re, "\\2", cl, perl = TRUE)
    from_aa  <- sub(change_re, "\\3", cl, perl = TRUE)  # e.g., "G 1.000" or "V 0.999"
    to_codon   <- sub(change_re, "\\4", cl, perl = TRUE)
    to_aa    <- sub(change_re, "\\5", cl, perl = TRUE)
    
    data.frame(
      branch  = b_id[i],
      species = species[i],
      coord   = coord,
      from_codon = from_codon,
      from_aa  = from_aa,
      to_codon   = to_codon,
      to_aa    = to_aa,
      change_line = cl,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, out)
}

library(Biostrings)

# Use built-in BLOSUM62 matrix
data(BLOSUM62)

BLOSUM62 = data.frame(BLOSUM62)



dirs <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)[-1]

for (d in dirs) {
  path <- file.path(d, "branches")
  
  if (file.exists(path) && file.info(path)$size > 0) {
    df <- parse_branch_changes(path)
    
    # Safety: ensure we have a data.frame with required cols; otherwise skip
    if (!is.data.frame(df)) next
    
    df <- df[, !(names(df) %in% c("branch", "from_codon", "to_codon", "change_line")), drop = FALSE]
    
    if (!all(c("from_aa","to_aa") %in% names(df)) || nrow(df) == 0) next
    
    df$from_aa <- sub(" .*", "", df$from_aa)
    df$to_aa   <- sub(" .*", "", df$to_aa)
    
    df <- df[df$from_aa != df$to_aa, , drop = FALSE]
    if (nrow(df) == 0) next
    
    df <- na.omit(df)
    if (nrow(df) == 0) next
    
    r <- match(df$from_aa, rownames(BLOSUM62))
    c <- match(df$to_aa,   colnames(BLOSUM62))
    df$score <- BLOSUM62[cbind(r, c)]
    
    out_path <- file.path(d, "branchchanges.csv")
    write.csv(df, out_path, row.names = FALSE)
  }
}

