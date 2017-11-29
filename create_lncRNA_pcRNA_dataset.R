library(stringr)
library(GeneStructureTools)
library(BSgenome.Hsapiens.UCSC.hg38)

getOrfs <- function(transcripts){
  transcripts <- transcripts[transcripts$type=="exon"]
  transcripts$exon_number <-
    as.numeric(transcripts$exon_number)
  order <-
    order(transcripts$transcript_id, transcripts$exon_number)
  transcripts <- transcripts[order]
  #transcripts$seq <-
  #  as.character(Biostrings::getSeq(g, transcripts))
  
  seqs_comb <-
    aggregate(seq ~ transcript_id, mcols(transcripts), function(x)
      (paste(x, collapse = "")))
  ids <- as.character(seqs_comb$transcript_id)
  seq_cat <- seqs_comb$seq
  rm <- which(grepl("N", seq_cat))
  
  if (length(rm) > 0) {
    seq_cat <- seq_cat[-rm]
    remove_id <- ids[rm]
    ids <- ids[-rm]
    transcripts <-
      transcripts[-which(transcripts$transcript_id %in% remove_id)]
  }
  
  # 3 frames
  seq_cat <-
    c(seq_cat, str_sub(seq_cat, 2), str_sub(seq_cat, 3))
  frames <- rep(c(1, 2, 3), each = length(ids))
  ids <- c(ids, ids, ids)
  
  orf <-
    suppressWarnings(unlist(lapply(seq_cat, function(x)
      as.character(Biostrings::translate(Biostrings::DNAString(x))))))
  
  orf_df <- data.frame(
    id = ids,
    aa_sequence = orf,
    frame = frames,
    stringsAsFactors = FALSE
  )
  
  orf_df$seq_length <- nchar(orf_df$aa_sequence)
  orf_df$seq_length_nt <- nchar(seq_cat) + orf_df$frame -1
  
  start_sites <- lapply(orf_df$aa_sequence, function(x) c(1, ((stringr::str_locate_all(x, "M")[[1]][,1]) -1 )))
  stop_sites <- lapply(orf_df$aa_sequence, function(x) c(((stringr::str_locate_all(x, "[*]")[[1]][,1])), nchar(x)))
  
  max_loc <-
    mapply(function(x, y)
      maxLocation(x, y), start_sites, stop_sites)
  
  orf_df$start_site <- max_loc[1, ]
  orf_df$stop_site <- max_loc[2, ]
  
  
  orf_df$orf_sequence <-
    stringr::str_sub(orf_df$aa_sequence, orf_df$start_site, orf_df$stop_site - 1)
  orf_df$orf_length <- nchar(orf_df$orf_sequence)
  
  orf_df$start_site_nt <-
    (orf_df$start_site * 3)- 3 + orf_df$frame
  orf_df$stop_site_nt <- (orf_df$orf_length * 3) + orf_df$start_site_nt + 3
  orf_df$utr3_length <-
    (orf_df$seq_length_nt - orf_df$stop_site_nt) + 1
  
  orf_df$aa_sequence <- NULL
  
  orf_df <- plyr::arrange(orf_df, plyr::desc(orf_length))
  orf_df <- orf_df[!duplicated(orf_df$id), ]
  
  orf_df <- plyr::arrange(orf_df, id)
  m <- match(orf_df$id, transcripts$transcript_id)
  orf_df$gene_id <- transcripts$gene_id[m]
  orf_df <- orf_df[,c(1, ncol(orf_df), 2:(ncol(orf_df)-1))]
  
  m <- match(orf_df$id, seqs_comb$transcript_id)
  orf_df$seq <- seqs_comb$seq[m]
  
  seq <- orf_df$seq
  sequence_n=lapply(seq, function(x) str_split(x,"")[[1]])
  
  transcript_lengths=nchar(seqs)
  
  sequence_pos_bias <- data.frame(A1=rep(NA, length(seq)),
                                  C1=rep(NA, length(seq)),
                                  G1=rep(NA, length(seq)),
                                  T1=rep(NA, length(seq)),
                                  A2=rep(NA, length(seq)),
                                  C2=rep(NA, length(seq)),
                                  G2=rep(NA, length(seq)),
                                  T2=rep(NA, length(seq)),
                                  A3=rep(NA, length(seq)),
                                  C3=rep(NA, length(seq)),
                                  G3=rep(NA, length(seq)),
                                  T3=rep(NA, length(seq)))
  
  pos1=seq(1, max(transcript_lengths),3)
  sequence_pos_bias$A1 <- unlist(lapply(sequence_n, function(x) sum(str_count(x[pos1], pattern="A"), na.rm=TRUE)))
  sequence_pos_bias$C1 <- unlist(lapply(sequence_n, function(x) sum(str_count(x[pos1], pattern="C"), na.rm=TRUE)))
  sequence_pos_bias$G1 <- unlist(lapply(sequence_n, function(x) sum(str_count(x[pos1], pattern="G"), na.rm=TRUE)))
  sequence_pos_bias$T1 <- unlist(lapply(sequence_n, function(x) sum(str_count(x[pos1], pattern="T"), na.rm=TRUE)))
  
  pos2=seq(2, max(transcript_lengths),3)
  sequence_pos_bias$A2 <- unlist(lapply(sequence_n, function(x) sum(str_count(x[pos2], pattern="A"), na.rm=TRUE)))
  sequence_pos_bias$C2 <- unlist(lapply(sequence_n, function(x) sum(str_count(x[pos2], pattern="C"), na.rm=TRUE)))
  sequence_pos_bias$G2 <- unlist(lapply(sequence_n, function(x) sum(str_count(x[pos2], pattern="G"), na.rm=TRUE)))
  sequence_pos_bias$T2 <- unlist(lapply(sequence_n, function(x) sum(str_count(x[pos2], pattern="T"), na.rm=TRUE)))
  
  pos3=seq(3, max(transcript_lengths),3)
  sequence_pos_bias$A3 <- unlist(lapply(sequence_n, function(x) sum(str_count(x[pos3], pattern="A"), na.rm=TRUE)))
  sequence_pos_bias$C3 <- unlist(lapply(sequence_n, function(x) sum(str_count(x[pos3], pattern="C"), na.rm=TRUE)))
  sequence_pos_bias$G3 <- unlist(lapply(sequence_n, function(x) sum(str_count(x[pos3], pattern="G"), na.rm=TRUE)))
  sequence_pos_bias$T3 <- unlist(lapply(sequence_n, function(x) sum(str_count(x[pos3], pattern="T"), na.rm=TRUE)))
  
  APOS=apply(sequence_pos_bias[,grep("A", colnames(sequence_pos_bias))], 1, function(x) (max(x)/(min(x)+1)))
  CPOS=apply(sequence_pos_bias[,grep("C", colnames(sequence_pos_bias))], 1, function(x) (max(x)/(min(x)+1)))
  GPOS=apply(sequence_pos_bias[,grep("G", colnames(sequence_pos_bias))], 1, function(x) (max(x)/(min(x)+1)))
  TPOS=apply(sequence_pos_bias[,grep("T", colnames(sequence_pos_bias))], 1, function(x) (max(x)/(min(x)+1)))
  
  orf_df <- cbind(orf_df, APOS, CPOS, GPOS, TPOS)
  
  
  number_C <- str_count(seq, c("C"))
  number_G <- str_count(seq, c("G"))
  number_T <- str_count(seq, c("T"))
  number_A <- str_count(seq, c("A"))
  
  GC_percent <- (number_C + number_G)/(number_A + number_C+ number_G +number_T)
  orf_df$GC_percent <- GC_percent
  orf_df$A_count <- number_A
  orf_df$C_count <- number_C
  orf_df$G_count <- number_G
  orf_df$T_count <- number_T
  
  orf_df$first_nt <- str_sub(seq, 1,1)
  orf_df$random_X <- sample(1:10000, nrow(orf_df), replace = TRUE)
  orf_df$random_Y <- rnorm(mean=10, sd=2, nrow(orf_df))
  
  return(orf_df)
}

# read in GTF
gencode <- rtracklayer::import("gencode.v27.annotation.gtf")

# add some extra annotations
gencode <- GeneStructureTools::addBroadTypes(gencode)
gencode <- GeneStructureTools::UTR2UTR53(gencode)

# extract "high quality" lncRNAs and protein coding transcripts

pcRNA_transcripts <- gencode[which(gencode$transcript_type=="protein_coding" & gencode$transcript_support_level < 2)]
lncRNA_transcripts <- gencode[which(gencode$transcript_type_broad=="lncRNA" & gencode$transcript_support_level < 2)]

# get sequences
g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
seqs <- Biostrings::getSeq(g, pcRNA_transcripts)
pcRNA_transcripts$seq <- as.character(seqs)
seqs <- Biostrings::getSeq(g, lncRNA_transcripts)
lncRNA_transcripts$seq <- as.character(seqs)

# get some sequence features (ORF and fickett scores)
orf_lncRNA <- getOrfs(lncRNA_transcripts)
orf_pcRNA <- getOrfs(pcRNA_transcripts)

orf_pcRNA$set <- "pcRNA"
orf_lncRNA$set <- "lncRNA"

rm_col <- match(c('orf_sequence','id','seq','gene_id'), colnames(orf_pcRNA))
all_data <- rbind(orf_pcRNA[,-rm_col], orf_lncRNA[,-rm_col])

all_data <- all_data[,c(1,3,2,7,4,8,5,6,9:ncol(all_data))]

write.csv(all_data, file='lncRNA_v_pcRNA.csv', row.names = FALSE)