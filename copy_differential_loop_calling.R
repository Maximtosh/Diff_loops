## Load packages
# if (!requireNamespace("remotes", quietly = TRUE))
#   install.packages("remotes")
# remotes::install_github("EricSDavis/hictoolsr")
library(hictoolsr)
library(dbscan)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("apeglm")
library("apeglm")


# p11_loops <- "/tank/projects/gmaxim_mouse/dots_8more_p11_new.tsv"
# p3_loops <- "/tank/projects/gmaxim_mouse/dots_8more_p3_new.tsv"
# w11_loops <- "/tank/projects/gmaxim_mouse/dots_8more_w11_new.tsv"
# w3_loops <- "/tank/projects/gmaxim_mouse/dots_8more_w3_new.tsv"

p11_loops <- "/tank/projects/gmaxim_mouse/p11_8d_loops.tsv"
p3_loops <- "/tank/projects/gmaxim_mouse/p3_8d_loops.tsv"
w11_loops <- "/tank/projects/gmaxim_mouse/w11_8d_loops.tsv"
w3_loops <- "/tank/projects/gmaxim_mouse/w3_8d_loops.tsv"

loops <-
  mergeBedpe(bedpeFiles = c(p11_loops, w11_loops, p3_loops, w3_loops), res = 25e3) |>
  as_ginteractions()

hicFiles <- 
  c("/tank/projects/gmaxim_mouse/PS19_11m.hic", 
    "/tank/projects/gmaxim_mouse/WT_11m.hic", 
    "/tank/projects/gmaxim_mouse/PS19_3m.hic", 
    "/tank/projects/gmaxim_mouse/WT_3m.hic")
loopCounts <- extractCounts(bedpe = loops,
                            hic = hicFiles,
                            chroms = c(paste0("chr", 1:19), "chrX"),
                            res = 25e3,
                            norm = 'NONE',
                            matrix = 'observed')
library(InteractionSet)
## Simplify column names
colnames(mcols(loopCounts)) <- 
  gsub(pattern = "(PS19|WT)_(11m|3m)\\.hic",
       replacement = "\\1_\\2",
       x = colnames(mcols(loopCounts)))
library(DESeq2)

## Isolate count matrix
cnts <- 
  mcols(loopCounts)[grep("WT|P", colnames(mcols(loopCounts)))] |>
  as.matrix()

## Create colData from column names
colData <- 
  do.call(rbind, strsplit(x = colnames(cnts), split = "_")) |>
  as.data.frame(stringsAsFactors = TRUE) |>
  `colnames<-`(value = c("condition", "age"))

## Build DESeq data set
dds <- 
  DESeqDataSetFromMatrix(countData = cnts,
                         colData = colData,
                         design = ~age + condition)
print(assay(dds,  normalized = TRUE))

res <- DESeq(dds)
print(res)
results(res, contrast =c('condition', 'WT', 'PS19'))
resu_ifc <- lfcShrink(res, coef = "condition_WT_vs_PS19", type="apeglm")
plotMA(resu_ifc)
print(resu_ifc)

plotPCA(vst(dds), intgroup = "condition") + ggplot2::theme(aspect.ratio = 1)

## Attach DESeq2 results
mcols(loopCounts) <- cbind(mcols(loopCounts), res)


### Normal loops (without 8 added anchors) ###

p11_loops <- "/tank/projects/gmaxim_mouse/test.dots_res_p11.25000.tsv"
p3_loops <- "/tank/projects/gmaxim_mouse/test.dots_res_p3.25000.tsv"
w11_loops <- "/tank/projects/gmaxim_mouse/test.dots_res_w11.25000.tsv"
w3_loops <- "/tank/projects/gmaxim_mouse/test.dots_res_w3.25000.tsv"

# p11_loops <- "/tank/projects/gmaxim_mouse/p11_25kb_nochr.csv"
# p3_loops <- "/tank/projects/gmaxim_mouse/p3_25kb_nochr.csv"
# w11_loops <- "/tank/projects/gmaxim_mouse/w11_25kb_nochr.csv"
# w3_loops <- "/tank/projects/gmaxim_mouse/w3_25kb_nochr.csv"

## Merge loops and convert to GInteractions
loops <-
  mergeBedpe(bedpeFiles = c(p11_loops, w11_loops, p3_loops, w3_loops), res = 25e3) |>
  as_ginteractions()

## Hi-C file paths from GEO
hicFiles <- 
  c("/tank/projects/gmaxim_mouse/PS19_11m.hic", 
    "/tank/projects/gmaxim_mouse/WT_11m.hic", 
    "/tank/projects/gmaxim_mouse/PS19_3m.hic", 
    "/tank/projects/gmaxim_mouse/WT_3m.hic")

## Extract Hi-C counts between loop pixels
loopCounts <- extractCounts(bedpe = loops,
                            hic = hicFiles,
                            chroms = c(paste0("chr", 1:19), "chrX"),
                            res = 25e3,
                            norm = 'NONE',
                            matrix = 'observed')
head(loops)
print(loopCounts)
# fwrite(loopCounts, "/tank/projects/gmaxim_mouse/loopcounts.tsv", sep = "\t", col.names = TRUE)


## Load package
library(InteractionSet)

## Simplify column names
colnames(mcols(loopCounts)) <- 
  gsub(pattern = "(PS19|WT)_(11m|3m)\\.hic",
       replacement = "\\1_\\2",
       x = colnames(mcols(loopCounts)))


head(loopCounts)

## Load package
library(DESeq2)

## Isolate count matrix
cnts <- 
  mcols(loopCounts)[grep("WT|P", colnames(mcols(loopCounts)))] |>
  as.matrix()

head(cnts)

## Create colData from column names
colData <- 
  do.call(rbind, strsplit(x = colnames(cnts), split = "_")) |>
  as.data.frame(stringsAsFactors = TRUE) |>
  `colnames<-`(value = c("condition", "age"))

colData

## Build DESeq data set
dds <- 
  DESeqDataSetFromMatrix(countData = cnts,
                         colData = colData,
                         design = ~age + condition)
print(assay(dds,  normalized = TRUE))
# a <- DESeq(dds)
# resultsNames(a)
## Run DEseq analysis
# res <-
#   DESeq(dds) |>
#   lfcShrink(coef = "condition_WT_vs_PS19", type="apeglm")

res <- DESeq(dds)
print(res)
results(res, contrast =c('condition', 'WT', 'PS19'))
resu_ifc <- lfcShrink(res, coef = "condition_WT_vs_PS19", type="apeglm")

summary(resu_ifc)  

plotMA(resu_ifc)
plotPCA(vst(dds), intgroup = "condition") + ggplot2::theme(aspect.ratio = 1)

r <- as.data.frame(resu_ifc)


loopCounts@regions

ranges1 <- anchors(loopCounts, type = "first")
ranges2 <- anchors(loopCounts, type = "second")

first_df <- as.data.frame(ranges1)
second_df <- as.data.frame(ranges2)

# Extract metadata columns as data.frame
meta_df <- as.data.frame(mcols(loopCounts))

combined_df <- data.frame(
  seqnames1 = first_df$seqnames,
  start1 = first_df$start,
  end1 = first_df$end,
  seqnames2 = second_df$seqnames,
  start2 = second_df$start,
  end2 = second_df$end,
  meta_df
)
head(combined_df[1:7])
fwrite(combined_df[1:7], "/tank/projects/gmaxim_mouse/combined_ranges.tsv", sep = "\t", col.names = TRUE)


## Attach DESeq2 results
mcols(loopCounts) <- cbind(mcols(loopCounts), res)

## Separate WT/FS-specific loops
wtLoops <- loopCounts[loopCounts$padj <= 0.01 &
                        loopCounts$log2FoldChange > 0]

pLoops <- loopCounts[loopCounts$padj <= 0.01 &
                        loopCounts$log2FoldChange < 0]

summary(wtLoops)
summary(pLoops)

## Hi-C file paths from GEO
# hicFiles <- c("bal_Hi-C_TauP301S_11months_25k.cool", "bal_Hi-C_TauP301S_3months_25k.cool")

## Load package
## Attach DESeq2 results
mcols(loopCounts) <- cbind(mcols(loopCounts), res)

## Separate WT/FS-specific loops
wtLoops <- loopCounts[loopCounts$padj <= 0.01 &
                         loopCounts$log2FoldChange > 0]

fsLoops <- loopCounts[loopCounts$padj <= 0.01 &
                         loopCounts$log2FoldChange < 0]

summary(wtLoops)
summary(fsLoops)




# EXAMPLE
## Load packages
library(hictoolsr)
library(dbscan)

## Define WT and FS loop file paths
wt_loops <- system.file("extdata/WT_5kbLoops.txt", package = "hictoolsr")
fs_loops <- system.file("extdata/FS_5kbLoops.txt", package = "hictoolsr")

## Merge loops and convert to GInteractions
loops_ex <- 
  mergeBedpe(bedpeFiles = c(wt_loops, fs_loops), res = 10e3) |>
  as_ginteractions()

head(loops)
## Load package
library(InteractionSet)

## Simplify column names
colnames(mcols(loopCounts)) <- 
  gsub(pattern = "GSM.*_IDR_(WT|FS)_A9_(1|2)_(1|2)_.*", 
       replacement = "\\1_\\2_\\3",
       x = colnames(mcols(loopCounts)))

head(loopCounts)
## Load package
library(DESeq2)

## Isolate count matrix
cnts_ex <- 
  mcols(loopCounts)[grep("WT|FS", colnames(mcols(loopCounts)))] |>
  as.matrix()

head(cnts_ex)
## Create colData from column names
colData_ex <- 
  do.call(rbind, strsplit(x = colnames(cnts_ex), split = "_")) |>
  as.data.frame(stringsAsFactors = TRUE) |>
  `colnames<-`(value = c("condition", "biorep", "techrep"))


## Build DESeq data set
dds_ex <- 
  DESeqDataSetFromMatrix(countData = cnts_ex,
                         colData = colData_ex,
                         design = ~techrep + biorep + condition)
dds_ex

## Run DEseq analysis
res_ex <-
  DESeq(dds_ex) |>
  lfcShrink(coef = "condition_WT_vs_FS", type="apeglm")

summary(res_ex)  
