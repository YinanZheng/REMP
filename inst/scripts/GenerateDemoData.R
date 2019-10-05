## Generate REMP demo data:

library(REMP)
library(devtools)

if(!isNamespaceLoaded("IlluminaHumanMethylation450kanno.ilmn12.hg19")) attachNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19")
if(!isNamespaceLoaded("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")) attachNamespace("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

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
ILMN450kEPIC.GR_38 = REMP::.liftOver_Hg19toHg38(ILMN450kEPIC.GR)

### Alu demo
Alu.GR <- fetchRMSK(REtype = "Alu", genome = "hg19", verbose = TRUE)
Alu38.GR <- fetchRMSK(REtype = "Alu", genome = "hg38", verbose = TRUE)

# Keep Alu only with ILMN profiled CpG
Alu_ILMN.GR = subsetByOverlaps(Alu.GR, ILMN450kEPIC.GR, ignore.strand = TRUE)
Alu38_ILMN.GR = subsetByOverlaps(Alu38.GR, ILMN450kEPIC.GR_38, ignore.strand = TRUE)

# 
RE_CpG <- findRECpG(RE = Alu_ILMN.GR, REtype = "Alu", genome = "hg19", verbose = TRUE)
RE38_CpG <- findRECpG(RE = Alu38_ILMN.GR, REtype = "Alu", genome = "hg38", verbose = TRUE)

RE_CpG_flanking <- .twoWayFlank(granges(RE_CpG), 900)
RE38_CpG_flanking <- .twoWayFlank(granges(RE38_CpG), 900)

HITS <- findOverlaps(RE_CpG_flanking, ILMN450kEPIC.GR, ignore.strand = TRUE)
HITS_38 <- findOverlaps(RE38_CpG_flanking, ILMN450kEPIC.GR_38, ignore.strand = TRUE)

## Part 1: RE CpG
RE_NeibCpG <- RE_CpG[queryHits(HITS), ]
RE38_NeibCpG <- RE38_CpG[queryHits(HITS_38), ]

## Total RE that can be predicted (contains neighboring CpGs within
## given window)

## Part 2: Neighboring ILMN CpG
RE_NeibCpG_ILMN <- ILMN450kEPIC.GR[subjectHits(HITS), ]
RE_NeibCpG_ILMN <- DataFrame(RE_NeibCpG_ILMN.GR = granges(RE_NeibCpG_ILMN), 
                             mcols(RE_NeibCpG_ILMN))
mcols(RE_NeibCpG) <- DataFrame(mcols(RE_NeibCpG), RE_NeibCpG_ILMN)

#
RE38_NeibCpG_ILMN <- ILMN450kEPIC.GR_38[subjectHits(HITS_38), ]
RE38_NeibCpG_ILMN <- DataFrame(RE_NeibCpG_ILMN.GR = granges(RE38_NeibCpG_ILMN), 
                             mcols(RE38_NeibCpG_ILMN))
mcols(RE38_NeibCpG) <- DataFrame(mcols(RE38_NeibCpG), RE38_NeibCpG_ILMN)

## Remove singleton (RE CpGs with only one neighboring ILMN CpG)
RE_NeibCpG$distance <- abs(start(RE_NeibCpG$RE_NeibCpG_ILMN.GR) - start(RE_NeibCpG))
RE_NeibCpG <- RE_NeibCpG[RE_NeibCpG$distance > 1,]
RE_NeibCpG_GR <- granges(RE_NeibCpG)
RE_NeibCpG <- subsetByOverlaps(RE_NeibCpG, unique(RE_NeibCpG_GR[duplicated(RE_NeibCpG_GR)]))

Alu_ILMN.GR <- Alu_ILMN.GR[Alu_ILMN.GR$Index %in% unique(RE_NeibCpG$Index)]

#
RE38_NeibCpG$distance <- abs(start(RE38_NeibCpG$RE_NeibCpG_ILMN.GR) - start(RE38_NeibCpG))
RE38_NeibCpG <- RE38_NeibCpG[RE38_NeibCpG$distance > 1,]
RE38_NeibCpG_GR <- granges(RE38_NeibCpG)
RE38_NeibCpG <- subsetByOverlaps(RE38_NeibCpG, unique(RE38_NeibCpG_GR[duplicated(RE38_NeibCpG_GR)]))

Alu38_ILMN.GR <- Alu38_ILMN.GR[Alu38_ILMN.GR$Index %in% unique(RE38_NeibCpG$Index)]


## Find same Alu between two build:
Alu_ILMN.GR$id <- paste0(seqnames(Alu_ILMN.GR), strand(Alu_ILMN.GR), Alu_ILMN.GR$name, Alu_ILMN.GR$score)
Alu38_ILMN.GR$id <- paste0(seqnames(Alu38_ILMN.GR), strand(Alu38_ILMN.GR), Alu38_ILMN.GR$name, Alu38_ILMN.GR$score)

commonid <- intersect(Alu_ILMN.GR$id, Alu38_ILMN.GR$id )

Alu_ILMN.GR_common <- Alu_ILMN.GR[match(commonid, Alu_ILMN.GR$id)]
Alu38_ILMN.GR_common <- Alu38_ILMN.GR[match(commonid, Alu38_ILMN.GR$id)]

## 
Alu38_ILMN.GR_common_LO <- REMP::.liftOver_Hg38toHg19(Alu38_ILMN.GR_common)
Alu_ILMN.GR_common_true <- subsetByOverlaps(Alu_ILMN.GR_common, Alu38_ILMN.GR_common_LO, type = "equal")

Alu_ILMN.GR_common_LO <- REMP::.liftOver_Hg19toHg38(Alu_ILMN.GR_common)
Alu38_ILMN.GR_common_true <- subsetByOverlaps(Alu38_ILMN.GR_common, Alu_ILMN.GR_common_LO, type = "equal")

identical(Alu_ILMN.GR_common_true$id, Alu38_ILMN.GR_common_true$id)
identical(Alu_ILMN.GR_common_true$name, Alu38_ILMN.GR_common_true$name)
identical(Alu_ILMN.GR_common_true$score, Alu38_ILMN.GR_common_true$score)


#
set.seed(2017)
ind <- sample(seq_len(length(Alu_ILMN.GR_common_true)), 500)

Alu.hg19.demo = GRanges(seqnames = factor(as.character(seqnames(Alu_ILMN.GR_common_true))[ind], levels = paste0("chr", 1:22)),
                   IRanges(start = as.integer(start(Alu_ILMN.GR_common_true))[ind],
                           end = as.integer(end(Alu_ILMN.GR_common_true))[ind]),
                   strand = as.character(strand(Alu_ILMN.GR_common_true))[ind],
                   name = as.character(Alu_ILMN.GR_common_true$name)[ind],
                   score = as.integer(Alu_ILMN.GR_common_true$score)[ind],
                   Index = Rle(as.character(Alu_ILMN.GR_common_true$Index)[ind]))

Alu.hg38.demo = GRanges(seqnames = factor(as.character(seqnames(Alu38_ILMN.GR_common_true))[ind], levels = paste0("chr", 1:22)),
                        IRanges(start = as.integer(start(Alu38_ILMN.GR_common_true))[ind],
                                end = as.integer(end(Alu38_ILMN.GR_common_true))[ind]),
                        strand = as.character(strand(Alu38_ILMN.GR_common_true))[ind],
                        name = as.character(Alu38_ILMN.GR_common_true$name)[ind],
                        score = as.integer(Alu38_ILMN.GR_common_true$score)[ind],
                        Index = Rle(as.character(Alu38_ILMN.GR_common_true$Index)[ind]))

Alu.hg19.demo <- sort(Alu.hg19.demo)
Alu.hg38.demo <- sort(Alu.hg38.demo)

Alu.hg38.demo_LO <- REMP::.liftOver_Hg38toHg19(Alu.hg38.demo)
Alu.hg19.demo_LO <- REMP::.liftOver_Hg19toHg38(Alu.hg19.demo)

identical(granges(Alu.hg38.demo_LO), granges(Alu.hg19.demo))
identical(granges(Alu.hg19.demo_LO), granges(Alu.hg38.demo))

#####################################
Alu.hg19.demo
Alu.hg38.demo

seqnames(Alu.hg19.demo)
seqnames(Alu.hg38.demo)

## Index must be Rle
Alu.hg19.demo$Index
Alu.hg38.demo$Index

use_data(Alu.hg19.demo, overwrite = TRUE)
use_data(Alu.hg38.demo, overwrite = TRUE)

