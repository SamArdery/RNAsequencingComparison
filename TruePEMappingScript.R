#!/usr/bin/env Rscript

.libPaths("/projects/munger-lab/ArderyProject/R/x86_64-pc-linux-gnu-library/4.1/")

library(tidyverse)
library(sva)
library(qtl2)

#load DO185_mESC_emse-m4_29180419.RData and PEMappingData.RData
load("/projects/munger-lab/ArderyProject/src/RNAseq_comp/MappingScripts/DO185_mESC_emase_m4_20180419.RData")
rm(expr, matches, raw.expr, DO185_mESC_emase_m4.cmd, exprNotes)
load("/projects/munger-lab/ArderyProject/src/RNAseq_comp/DataCollection/PEMappingData.RData")
load("/projects/munger-lab/ArderyProject/src/RNAseq_comp/MappingScripts/DO_mESC_eQTL_peaks.RData")
rm(peaks, cmd, DO185_mESC_emase_m4.cmd)
n.cores <- 16
final <- PEfinal
PEDatacmd <- "This expression data was obtained by running the DataCollection.Rmd file"
outdir <- "/projects/munger-lab/ArderyProject/eQTLRemap"

split_map <- function (map, chr_names = NULL)
{
  map <- reorder_map_table(map, chr_names = chr_names)
  pos <- as.numeric(map[, 2])
  chr <- map[, 1]
  uchr <- unique(chr)
  names(pos) <- rownames(map)
  lapply(split(pos, factor(chr, uchr)), sort)
}

#This calls reorder_map_table() and the code for that is:
reorder_map_table <- function (map_tab, chr_col = 1, pos_col = 2, chr_names = NULL)
{
  chr <- map_tab[, chr_col]
  if (is.null(chr_names))
    chr_names <- unique(chr)
  chr <- factor(chr, levels = chr_names)
  pos <- suppressWarnings(as.numeric(map_tab[, pos_col]))
  map_tab[order(chr, pos, seq_len(nrow(map_tab))), , drop = FALSE]
}

#Taking the raw expression data and creating a normalized matrix
finalmat <- dplyr::select(final, -sample, -rep) %>%
  spread(id, count) %>%
  remove_rownames() %>%
  as.data.frame() %>%
  column_to_rownames(var="locus") %>%
  as.matrix()
q75 <- apply(finalmat, 2, quantile, 0.75)
ratio <- mean(q75)/q75
finalmat2 <- sweep(finalmat, 2, ratio, FUN="*")

raw.expr.paired <- finalmat
norm.expr.paired <- finalmat2

dat <- log(norm.expr.paired + 1)
mod <- model.matrix(~ sex, data=covarTidy)
exprComBat <- ComBat(dat=dat, batch=covarTidy$libraryprep, mod=mod, 
                     par.prior=TRUE, prior.plots=FALSE)
covar <- covar[,"sex", drop=FALSE]

# Prepare for mapping
rankZ <- function (x) {
  x <- rank(x, na.last = "keep", ties.method = "average")/(sum(!is.na(x)) + 1)
  qnorm(x)
}
exprZ <- apply(exprComBat, 1, rankZ)    # samples in rows, genes (==phenotype) in columns

# let's modify genoprobs
# This function is from Dan:
message("converting probs to qtl2 and calculating kinship matrix")
uchroms <- unique(map_dat$chr)
probs_3d_to_qtl2 <- function(probs) {
  # Convert to qtl2 genoprobs format
  # Similar to qtl2convert::probs_doqtl_to_qtl2()
  markers <- dimnames(probs)[[3]]
  chroms <- sapply(strsplit(markers, "_"), "[[", 1)
  newprobs <- vector("list", length(uchroms))
  names(newprobs) <- uchroms
  for (chrom in uchroms) newprobs[[chrom]] <- probs[, , chroms == chrom]
  attr(newprobs, "is_x_chr") <- c(rep(FALSE, length(uchroms)-1),TRUE)
  attr(newprobs, "crosstype") <- "DO"
  attr(newprobs, "alleles") <- c("A","B","C","D","E","F","G","H")
  attr(newprobs, "alleleprobs") <- TRUE
  class(newprobs) <- c("calc_genoprob", "list")
  newprobs
}

genoprobs <- probs_3d_to_qtl2(probs)
# probs <- qtl2convert::probs_doqtl_to_qtl2(probs, map_dat,
#   chr_column="chr", pos_column="pos", marker_column="marker")
# kinship_loco <- qtl2geno::calc_kinship(probs, "loco", cores=n.cores)
kinship_loco <- qtl2::calc_kinship(genoprobs, "loco", cores=n.cores)

# For saving the pmap as well:
map_dat2 <- map_dat %>%
  separate(marker, into=c('chrom', 'pos_bp'), convert=T, remove=F) %>%
  mutate(n=1:n()) %>% as_tibble()
pmap <- split_map(dplyr::select(map_dat2, marker,
                                chr, pos_bp) %>% as.data.frame() %>%
                    tibble::remove_rownames() %>%
                    tibble::column_to_rownames('marker'))

#saving information necessary to perform mapping in interactive sessions
save(map_dat2, pmap, gmap, exprZ, kinship_loco, exprComBat, genoprobs, covar, file = "/projects/munger-lab/ArderyProject/src/RNAseq_comp/MappingScripts/TruePEforMapping.RData")

message(date(), " doing qtl mapping")
batchmap <- function(nbatch, thrA=5, thrX=5, ...) {
  # Do mapping in batches. I will only save significant peaks, for efficiency.
  # I am saving the scan files for comparison plots later.
  nn <- ncol(exprZ)
  ss <- round(seq(0, nn, length.out=nbatch + 1))
  peaks <- list()
  for (i in 1:nbatch) {
    start <- ss[i] + 1
    end <- ss[i + 1]
    cat(sprintf("batch %d: %d-%d\n", i, start, end))
    out <- qtl2::scan1(genoprobs, exprZ[, start:end, drop=FALSE],
                       kinship_loco, addcovar=covar, cores=n.cores, ...)
    scan.file <-  paste0(outdir, "/DO_mESC_eQTL_scan_",i,".RData")
    save(out, file=scan.file)
    peaks[[i]] <- qtl2::find_peaks(out, gmap, drop=1.5,
                                   threshold=thrA, thresholdX=thrX)   # returns a long & tidy dataset
  }
  do.call('rbind', peaks) %>% dplyr::select(-lodindex) %>%
    dplyr::rename(phenotype=lodcolumn, peak_chr=chr, peak_cM=pos)
}
#exprZ <- exprZ[, 1:10]   # testing...
n.batches <- max(c(round(ncol(exprZ)/1000), 2))
message("Mapping ", ncol(exprZ), " gene expression levels. Running in ", n.batches, " batches")
peaks <- batchmap(n.batches)

cmd <- paste0("This mapping data was created using TruePEMappingScript.R on ", Sys.Date())

outfile <- paste0(outdir, "/DO_mESC_paired_eQTL_peaks.RData")
save(peaks, cmd, map_dat, gmap, PEDatacmd, file=outfile)
message(date(), " finished")

