#
# Copyright 2023 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@iit.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

args = commandArgs(trailingOnly=TRUE)

for(v in args)
{
  vTmp <- strsplit(v,"=")[[1]]
  assign(vTmp[[1]],vTmp[[2]])
}

#load Biostrings package
library("Biostrings")

Run_MSSPE = function(mfa_file, primers_file, ovlp_window_size, search_window_size, kmer_size, num_acc_miss, num_max_it) {
  #import the multiple sequence alignment as a DNAStringSet object and check all sequences in the gapped alignment have the same length
  mfa <- readDNAStringSet(filepath = mfa_file, format = "fasta")
  mfa_width <- unique(width(mfa))
  if (length(mfa_width) > 1) {
    stop("There were issues with the Multiple Sequence Alignment file!\n")
  }
  
  #Checks
  if (length(ovlp_window_size) > mfa_width) {
    stop("ovlp_window_size is greater than the length of sequences in fasta file. Please check ovlp_window_size parameter!\n")
  }
  if (length(ovlp_window_size) < 2*search_window_size) {
    stop("ovlp_window_size is smaller than 2*search_window_size. Please consider decreasing ovlp_window_size or increasing search_window_size parameters!\n")
  }
  
  #identify windows of width 2*ovlp_window_size with overlap of ovlp_window_size
  cuts <- unique(c(seq(from = 1, to = mfa_width, by = ovlp_window_size), mfa_width))
  cuts_start <- cuts[seq_along(cuts) %% 2 == 1]
  cuts_end <- cuts[seq_along(cuts) %% 2 == 0]
  cuts_GRanges <- IRanges(start = cuts_start, end = cuts_end)
  #identify the first and last search_window_size bases for primers search
  cuts_GRanges_start <- resize(cuts_GRanges, width = search_window_size, fix = "start")
  cuts_GRanges_end <- resize(cuts_GRanges, width = search_window_size, fix = "end")
  #extract subsequences corresponding to selected window
  mfa_split_start <- extractAt(mfa, cuts_GRanges_start)
  mfa_split_end <- extractAt(mfa, cuts_GRanges_end)
  #identify all possible kmers of length kmer_size
  kmers_start <- 1:(search_window_size - kmer_size + 1)
  kmers_end <- kmer_size:search_window_size
  kmers_GRanges <- IRanges(start = kmers_start, end = kmers_end)
  #extract subsequences corresponding to all possible kmers of length kmer_size from the first search_window_size of each window
  mfa_split_seq_tmp <- lapply(mfa_split_start, function(x) extractAt(x, kmers_GRanges))
  mfa_split_seq <- lapply(mfa_split_seq_tmp, function(x) as.character(unlist(x)))
  #extract subsequences corresponding to the reverse complement of all possible kmers of length kmer_size from the last search_window_size of each window
  mfa_split_seq_rc_tmp <- lapply(mfa_split_end, function(x) reverseComplement(unlist(extractAt(x, kmers_GRanges))))
  mfa_split_seq_rc <- lapply(mfa_split_seq_rc_tmp, function(x) as.character(x))
  names(mfa_split_seq_rc) <- paste0(names(mfa_split_seq_rc), "_revComp")
  #merge all kmers and remove unwanted patterns in the alignment
  kmers_tmp <- c(mfa_split_seq, mfa_split_seq_rc)
  kmers_split <- split(kmers_tmp, names(kmers_tmp))
  kmers <- lapply(kmers_split, function(x) {
    tmp <- unname(unlist(x))
    ind_noCons <- grep(x = tmp, pattern = "-|N|R|Y|S|W|K|M|B|D|H|V|N")
    if (length(ind_noCons) > 0) {
      tmp <- tmp[-ind_noCons]
    }
    return(tmp)
  })
  
  #reiterate on the list of segments until specific criteria are met
  #at each iteration identify the most frequent kmer and remove corresponding segment, together with its reverse complement
  residual_kmers <- kmers
  names(residual_kmers) <- gsub(x = names(residual_kmers), pattern = "\\(", replacement = "_")
  names(residual_kmers) <- gsub(x = names(residual_kmers), pattern = "\\)", replacement = "_")
  candidate_primers <- c()
  names_primer_found <- list()
  counter <- 1
  while (length(residual_kmers) > num_acc_miss & counter < num_max_it) {
    candidate_primer_curr <- names(sort(table(unlist(residual_kmers)), decreasing = TRUE)[1])
    candidate_primers <- c(candidate_primers, candidate_primer_curr)
    ind_primer_found <- which(grepl(pattern = candidate_primer_curr, residual_kmers))
    segments_primer_found <- which(grepl(x = names(residual_kmers), pattern = paste0(unique(gsub(x = names(residual_kmers)[ind_primer_found], pattern = "_revComp", replacement = "")), collapse = "|")))
    names_primer_found[[counter]] <- names(residual_kmers)[ind_primer_found]
    residual_kmers <- residual_kmers[-segments_primer_found]
    counter <- counter + 1
  }
  
  #assign names to primers based on the amount of sequences in the database sharing the primer
  primers_names <- unlist(lapply(names_primer_found, function(x) {
    num_rc <- length(grep(pattern = "_revComp", x = x))
    num_fw <- length(x) - num_rc
    if (num_rc > 0 && num_fw == 0) {
      name <- paste0("Primer_rv_", num_rc, "_sequences")
    } else if (num_rc == 0 && num_fw > 0) {
      name <- paste0("Primer_fw_", num_fw, "_sequences")
    } else {
      name <- paste0("Primer_fw_", num_fw, "_sequences_", num_rc, "_rv_sequences")  
    }
  }))
  
  candidate_primers <- DNAStringSet(candidate_primers)
  names(candidate_primers) <- primers_names
  
  #filter candidate primers for homopolymers
  pattern_hompolymers <- c("AAAAA", "TTTTT", "CCCCC", "GGGGG")
  ind_homopol <- c()
  for (p in pattern_hompolymers) {
    ind_homopol <- c(ind_homopol, which(unlist(lapply(vmatchPattern(pattern = p, subject = candidate_primers), function(x) length(x) > 0))))
  }
  ind_homopol <- unique(ind_homopol)
  
  #filter candidate primers for melting temperature Tm
  GC_count <- unlist(lapply(vmatchPattern(pattern = "S", subject = candidate_primers, fixed = "subject"), function(x) sum(width(x))))
  AT_count <- unlist(lapply(vmatchPattern(pattern = "W", subject = candidate_primers, fixed = "subject"), function(x) sum(width(x))))
  Tm <- 64.9 + 41*(GC_count - 16.4)/(GC_count + AT_count)
  Tm_threshold <- mean(Tm) + 2*sd(Tm)
  ind_high_Tm <- which(Tm > Tm_threshold)
  
  #discard primers that do not met the criteria
  ind_discard <- c(ind_high_Tm, ind_homopol)
  if (length(ind_discard) > 0) {
    primers <- candidate_primers[-ind_discard]  
  } else {
    primers <- candidate_primers
  }
  
  if (length(primers) > 0) {
    #save primers to file
    writeXStringSet(x = primers, filepath = primers_file)
  } else {
    cat(sprintf("No primers found!\n"))
  }
}

#Run MSSPE
Run_MSSPE(mfa_file, primers_file, ovlp_window_size, search_window_size, kmer_size, num_acc_miss, num_max_it)
