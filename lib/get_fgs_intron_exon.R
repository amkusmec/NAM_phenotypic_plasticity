###
# Original script by Jinliang Yang

###################
get_fgs_intron_exon <- function(FGSgff="~/gxe-gwas2/data/ZmB73_5b_FGS.gff",
                                FGSinfo="~/gxe-gwas2/data/ZmB73_5b_FGS_info.txt"){
  ### read in the FGS_v2 gff file
  fgs <- read.table(FGSgff)
  names(fgs) <- c("seqname", "source", "feature", "start", "end", "score",
                  "strand", "frame", "attribute")
  message(sprintf("[ %s ] lines of FGSv2 were loaded!", nrow(fgs)))
  
  ### get canonical transcripts!
  geneinfo <- read.table(FGSinfo, header=TRUE)
  canon <- subset(geneinfo, is_canonical=="yes")
  #[1] 39656    12
  message(sprintf("[ %s ] canonical trascripts were loaded!", nrow(canon)))
  
  gene <- subset(fgs, feature %in% "gene")
  gene$geneid <- gsub(";Name=.*", "", gene$attribute)
  gene$geneid <- gsub("ID=", "", gene$geneid)
  message(sprintf("[ %s ] genes' [ %s ]  canonical trascripts were found!", 
                  nrow(gene), length(unique(gene$geneid)) ))
  
  exon <- subset(fgs, feature %in% "exon")
  exon$txid <- gsub(";Name=.*", "", exon$attribute)
  exon$txid <- gsub("Parent=", "", exon$txid)
  exon <- subset(exon, txid %in% canon$transcript_id)
  message(sprintf("[ %s ] genes' [ %s ]  canonical exons were found!", 
                  length(unique(exon$txid)), nrow(exon) ))
  
  intron <- subset(fgs, feature %in% "intron")
  intron$txid <- gsub(";Name=.*", "", intron$attribute)
  intron$txid <- gsub("Parent=", "", intron$txid)
  intron <- subset(intron, txid %in% canon$transcript_id)
  message(sprintf("[ %s ] genes' [ %s ]  canonical introns were found!", 
                  length(unique(intron$txid)), nrow(intron) ))
  
  outls <- list();
  outls[['gene']] <- gene;
  outls[['exon']] <- exon;
  outls[['intron']] <- intron;
  message(sprintf("A list object was output!"))
  return(outls)
  
}
############
#fgs <- get_fgs_intron_exon()