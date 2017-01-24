## Generate REMP demo data:

library(REMP)
library(devtools)

annot450k = minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotEPIC = minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

## Remove "ch" probes
annot450k = subset(annot450k, substring(Name, 1,2) != "ch")
annotEPIC = subset(annotEPIC, substring(Name, 1,2) != "ch")

## Remove sex probes
annot450k = subset(annot450k, !chr %in% c("chrX","chrY"))
annotEPIC = subset(annotEPIC, !chr %in% c("chrX","chrY"))

ILMN450k.GR = GenomicRanges::GRanges(seqnames = annot450k$chr, 
                                     ranges = IRanges::IRanges(start = annot450k$pos, end = annot450k$pos), 
                                     strand = annot450k$strand, 
                                     ILMNID = rownames(annot450k))

ILMNEPIC.GR = GenomicRanges::GRanges(seqnames = annotEPIC$chr,
                                     ranges = IRanges::IRanges(start = annotEPIC$pos, end = annotEPIC$pos),
                                     strand = annotEPIC$strand,
                                     ILMNID = rownames(annotEPIC))

ILMN450kEPIC.GR = subsetByOverlaps(ILMN450k.GR, ILMNEPIC.GR, ignore.strand = TRUE)
### Alu.demo
ah = AnnotationHub::AnnotationHub()
Alu.GR <- fetchRMSK(ah, "Alu", verbose = TRUE)

# Keep Alu only with ILMN profiled CpG
Alu_ILMN.GR = subsetByOverlaps(Alu.GR, ILMN450kEPIC.GR, ignore.strand = TRUE)

# 
RE_CpG <- findRECpG(Alu_ILMN.GR, "Alu", verbose = TRUE)

RE_CpG_flanking <- .twoWayFlank(granges(RE_CpG), 900)
HITS <- findOverlaps(RE_CpG_flanking, ILMN450kEPIC.GR, ignore.strand = TRUE)

## Part 1: RE CpG
RE_NeibCpG <- RE_CpG[queryHits(HITS), ]
## Total RE that can be predicted (contains neighboring CpGs within
## given window)

## Part 2: Neighboring ILMN CpG
RE_NeibCpG_ILMN <- ILMN450kEPIC.GR[subjectHits(HITS), ]
RE_NeibCpG_ILMN <- DataFrame(RE_NeibCpG_ILMN.GR = granges(RE_NeibCpG_ILMN), 
                             mcols(RE_NeibCpG_ILMN))
mcols(RE_NeibCpG) <- DataFrame(mcols(RE_NeibCpG), RE_NeibCpG_ILMN)

## Remove singleton (RE CpGs with only one neighboring ILMN CpG)
RE_NeibCpG$distance <- abs(start(RE_NeibCpG$RE_NeibCpG_ILMN.GR) - start(RE_NeibCpG))
RE_NeibCpG <- RE_NeibCpG[RE_NeibCpG$distance > 1,]
RE_NeibCpG_GR <- granges(RE_NeibCpG)
RE_NeibCpG <- subsetByOverlaps(RE_NeibCpG, unique(RE_NeibCpG_GR[duplicated(RE_NeibCpG_GR)]))

Alu_ILMN.GR <- Alu_ILMN.GR[Alu_ILMN.GR$Index %in% unique(RE_NeibCpG$Index)]

set.seed(2017)
Alu.demo = Alu_ILMN.GR[sample(seq_len(length(Alu_ILMN.GR)), 500),]
nrun(Alu.demo$Index)

use_data(Alu.demo, overwrite = TRUE)

## REMPset.demo
# initREMP(arrayType = "450k", REtype = "Alu", RE = Alu.demo, ncore = 1, verbose = TRUE)
# initREMP(arrayType = "EPIC", REtype = "Alu", RE = Alu.demo, ncore = 1, verbose = TRUE)
# 
# GM12878_450k <- getGM12878("450k") # Get GM12878 methylation data (450k array)
# GM12878_EPIC <- getGM12878("EPIC") # Get GM12878 methylation data (450k array)
# 
# REMPset.450k.demo <- remp(GM12878_450k, REtype = "Alu", autoTune = FALSE, param = 6, ncore = 1, verbose = TRUE)
# REMPset.EPIC.demo <- remp(GM12878_EPIC, REtype = "Alu", autoTune = FALSE, param = 6, ncore = 1, verbose = TRUE)
# 
# REMPset.450k.demo
# REMPset.EPIC.demo

# use_data(REMPset.450k.demo, overwrite = TRUE)


