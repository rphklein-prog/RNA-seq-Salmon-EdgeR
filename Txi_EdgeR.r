#!/usr/bin/env Rscript

# arguments: first is transcriptome name (ex EnsDb.Hsapiens.v75), the second is a file of 
# the data groups (example: "Control", "Control", "Treatment", "Treatment"), the third 
# through second to last are salmon quant file names, last is the name for the output file

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 4) {
stop("Usage: Rscript script.R <transcriptome_name> <data_groups_file> <salmon_files...> <output_file>")
}

transcriptome_name <- args[1]
DataGroups <- scan(args[2], what = "", sep = ",", quote = "\"", strip.white = TRUE)
files <- args[3:(length(args)-1)]
output_file <- args[length(args)]

#import all necessary libraries

suppressPackageStartupMessages({
library(tximportData)
library(ensembldb)
library(transcriptome_name, character.only = TRUE)
library(tximport)
library(edgeR)
})


# get names of samples and create table with a column for each sample's quant data for each
# transcript the transcript file
names(files) <- paste0("sample", seq_along(files))
edb <- get(transcriptome_name)
Tx <- transcripts(edb, return.type="DataFrame")
tx2gene <- subset(Tx, select=c("tx_name", "gene_id"))
txi.salmon <- tximport(files, type = "salmon", ignoreTxVersion=TRUE, tx2gene = tx2gene)

# create a DGEList object
cts <- txi.salmon$counts
x <- DGEList(cts,group=factor(DataGroups))
x.full <- x

# remove data with low counts in 2 or more samples
print(apply(x$counts, 2, sum))
keep <- rowSums(cpm(x)>10) >=2
x <- x[keep,]
x$samples$lib.size <- colSums(x$counts)

#calculate normalization factors
x <-calcNormFactors(x)
x <- estimateDisp(x)

#run glm with quasi-likelihood
x1 <- glmQLFit(x)
x2 <- glmQLFTest(x1)

#return top tags and write datatable
top <- topTags(x2, n=Inf)$table
write.table(top, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
