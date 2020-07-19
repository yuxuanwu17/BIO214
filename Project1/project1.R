library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
genome <- BSgenome.Hsapiens.UCSC.hg19
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

transcripts <- exonsBy(txdb, by="tx", use.names =TRUE)
# seqlevels(transcripts)
transcripts <- transcripts[seqnames(transcripts) == "chr21"]
tx_width <- sum(width(transcripts))
transcripts <- transcripts[tx_width>100 & tx_width<=500]

tx_seqs <- extractTranscriptSeqs(Hsapiens,transcripts)

writeXStringSet(tx_seqs,"tx_seqs.fasta")

system("RNAfold --temp=70 tx_seqs.fasta > RNAstructure.txt")
RNAfold_out <- readLines("RNAstructure.txt")
Struc_pred <- RNAfold_out[seq(3, by = 3, length.out = length(tx_seqs))]

## Remove the energy scores attached at the end:
Struc_pred <- gsub(" .*", "", Struc_pred)
Struc_pred <- BStringSet(Struc_pred)
names(Struc_pred) <- names(tx_seqs)
Struc_pred



## Construct the GRanges object for hybridized and nonhybridized regions:

## === Hint code, fill the xs  ============= ##
#Hyb_irl <- lapply(vmatchPattern("xxx",Struc_pred), reduce)
#nonHyb_irl <- lapply(Hyb_irl, xxxx) # Note: the `xxxx` is one of the inter-range method, please see BIO214_refcard.rmd.
## ===== Enter your code below  ============ ##

nonHyb_irl  <- lapply(vmatchPattern(".",Struc_pred), reduce)
Hyb_irl <- lapply(nonHyb_irl, gaps)

##===== Your Code is finished till here =====##

##Convert the IrangesList into GRanges

irl2grl <- function(irl) GRangesList( mapply(function(x,y) GRanges(seqnames = y,
                                                                   ranges = x,
                                                                   strand = "*"),irl,names(irl)) )

Hyb_gr <- unlist(irl2grl(Hyb_irl))
nonHyb_gr <- unlist(irl2grl(nonHyb_irl))

Hyb_gr ##The Granges for Hybridized regions on the transcript
nonHyb_gr ##The Granges for non-Hybridized (Looped) regions on the transcript

seqnames(Hyb_gr)
seqnames(nonHyb_gr)

Hyb_gr <- mapFromTranscripts(Hyb_gr,transcripts)
nonHyb_gr <- mapFromTranscripts(nonHyb_gr,transcripts)
##===== Your Code is finished till here =====##
seqnames(Hyb_gr)
seqnames(nonHyb_gr)

Hyb_GC <- letterFrequency(DNAStringSet(Views(Hyb_gr,start = )), c("G,C"), as.prob=TRUE)

Views(genome, Hyb_gr)

Views(Hyb_gr)
