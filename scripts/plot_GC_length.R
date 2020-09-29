#!/usr/bin/Rscript
# PAUVERT Charlie
# 2018-02-07
# Distribution of GC and length from a fasta file


# Fetch command line argument
args <- commandArgs(trailingOnly = TRUE)

# Read FASTA sequence
library(seqinr)
x<-read.fasta(args[1])

# Get data
library(plyr)
dat<-data.frame(
	Length = ldply(x,function(seq) summary(seq)$length)$V1,
	GC = ldply(x,function(seq) summary(seq)$GC)$V1
)


# Plot
library(ggplot2)
library(cowplot)

pl<-ggplot(dat, aes(x = Length,y = GC)) +
 geom_point(shape = 1) +
 geom_rug(color = "forestgreen",alpha=.5) +
 labs( title = paste(nrow(dat),"sequences"),
       x = "Length (bp)",y="G+C percentage") + 
 theme_cowplot() + 
 theme(text = element_text(face="bold"))


# Save locally
message(paste("File",paste0(args[1],"_length_GC.png"), "saved"))
ggsave(pl, file = paste0(args[1],"_length_GC.png"),width=4,height=4,dpi=300)
