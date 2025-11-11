#!/usr/bin/env Rscript

# this script takes RNA-seq data analyzed with EdgeR and uses GOseq to run Gene ontology 
# enrichment analysis. Assumes 2 conditions are being tested (ex- treatment and control). 
# Script takes 4 arguments: first is the topTags file, the second is the genome build 
# (ex- hg19), the third is the annotation database (ex- ensGene), the last is the output 
# file name

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 4) {
stop("Usage: Rscript script.R <edgeR_data> <genome_build> <transcriptome_name> <output_name")
}

## start with file from EdgeR
et12 <- read.table(args[1], header = TRUE, row.names = 1)
genome <- args[2]
database <- args[3]
output_name <- args[4]

## install the goseq program and annotation database file in R, you will only need to do 
# this part once

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) 
    BiocManager::install("org.Hs.eg.db")

if (!requireNamespace("goseq", quietly = TRUE)) 
	BiocManager::install("goseq")

## load the libraries
library(goseq)
library(org.Hs.eg.db)
library(GO.db)

## choose the significantly differentially expressed genes (up and down)and put them in a 
## new table

genes=as.integer(p.adjust(et12$table$PValue[et12$table$logFC!=0],method="BH")<.05)
names(genes)=row.names(et12$table[et12$table$logFC!=0,])

## identify the species in the samples
supportedOrganisms()[supportedOrganisms()$Genome==genome,]
pwf=nullp(genes,genome,database)
GO.wall=goseq(pwf,genome,database)

## correct for multiple hypothesis testing
enriched.GO <- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "BH") < .05]

## write the data to a table and save it
write.table(GO.wall, file=output_name, "_all.txt", sep="\t")

## repeat the process, but choose the significantly downregulated genes and put them in a 
## new table

genes=as.integer(p.adjust(et12$table$PValue[et12$table$logFC<0],method="BH")<.05)
names(genes)=row.names(et12$table[et12$table$logFC<0,])

## identify the species in the samples
> supportedOrganisms()[supportedOrganisms()$Genome==genome,]
pwf=nullp(genes,genome,database)
GO.wall=goseq(pwf,genome,database)

## correct for multiple hypothesis testing
enriched.GO <- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "BH") < .05]

## write the data to a table and save it
> write.table(GO.wall, file=output_name, "_downreg.txt", sep="\t")

## this time choose the significantly upregulated genes

genes=as.integer(p.adjust(et12$table$PValue[et12$table$logFC>0],method="BH")<.05)
names(genes)=row.names(et12$table[et12$table$logFC>0,])

## identify the species in the samples
> supportedOrganisms()[supportedOrganisms()$Genome==genome,]
pwf=nullp(genes,genome,database)
GO.wall=goseq(pwf,genome,database)

## correct for multiple hypothesis testing
enriched.GO <- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "BH") < .05]

## write the data to a table and save it
> write.table(GO.wall, file=output_name, "_upreg.txt", sep="\t")

