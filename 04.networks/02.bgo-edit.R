### Process BINGO output files from Cytoscape. This part of the analysis will
### have to be done by hand.

source("../00.load-packages.R")

for (o in c("bp", "mf", "cc")) {
  files <- list.files(paste0("bgo/", o), "*.bgo", full.names = TRUE)
  for (f in seq_along(files)) {
    go <- read_lines(files[f])
    go <- go[grep("GO-ID", go):length(go)]
    go <- str_split(go, "\t")
    
    temp <- matrix("", ncol = 9)
    if (length(go) > 1) {
      for (j in 2:length(go)) {
        temp <- rbind(temp, go[[j]])
      }
      ### This nasty expression is necessary because temp[-1, ] will return a
      ### vector if nrow(temp) == 2. A vector passed to as_tibble() is turned
      ### into a single column tbl_df.
      temp <- as_tibble(matrix(temp[-1, ], ncol = 9,
                               dimnames = list(NULL, make.names(go[[1]])))) %>%
        mutate(GO.ID = as.integer(GO.ID),
               p.value = as.numeric(p.value),
               corr.p.value = as.numeric(corr.p.value),
               x = as.integer(x),
               n = as.integer(n),
               X = as.integer(X),
               N = as.integer(N))
      
      stem <- unlist(str_split(files[f], "/"))[3]
      write_tsv(temp, path = paste0("bgo-edit/", o, "/", stem))
    }
  }
}
