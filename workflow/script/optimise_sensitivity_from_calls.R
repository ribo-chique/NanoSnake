#!/usr/bin/env Rscript

# This script carries out iterative sensitivity/specificity measurements from nanopore modification probability distribution data of IVT negative control RNA, using modkit sample-probs outputs probabilities.tsv and thresholds.tsv
# The user defines a desired sensitivity score, then this script finds the lowest probability threshold which satisfy this score
# Using thresholds.tsv, we can obtain the default 10% quantile probability score threshold
# Then by incrementally changing the thresholds to a higher value, the number of "passed" nucleotides will change
# The "accuracy" score is based on the number of passed unmodified bases divided by total number of passed bases. This is the same method as modkit validate.
# Modkit validate uses the entire dataset to compute sensitivity, while this relies on sample-probs which uses fewer sampled reads, but a much quicker execution with the ability to iterate.

# supress warnings
options(warn=-1)

suppressPackageStartupMessages(library(tidyverse))

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(parallel))

# option parsing #

parser <- OptionParser()

parser <- add_option(parser, 
                     opt_str = c("-i", "--input"), 
                     type = "character",
                     dest = 'input_dir',
                     help="Input directory (required). A full path to a directory containing subdirectories with a modkit sample-probs tsv. The parent directory is searched recursively."
)

parser <- add_option(parser, 
                     opt_str = c("-o", "--output"), 
                     type = "character",
                     dest = 'output_dir',
                     default=NA,
                     help="Output directory (optional, default is the same as the input directory). The output directory will contain sub-directories that contain variables associated with each input tsv file."
)

parser <- add_option(parser, 
                     opt_str = c("-a", "--accuracy"), 
                     type = "double",
                     dest = 'target_accuracy',
                     help="Desired accuracy score in percentages, which the output filter threshold will satisfy (required). eg 99 for 99% accuracy."
)

parser <- add_option(parser, 
                     opt_str = c("-p", "--processors"), 
                     type="integer", 
                     default=1,
                     dest = 'cores',
                     help="Number of processors to use for the anlaysis (default 1). Only helps when you are analysing more than one file at a time"
)

opt = parse_args(parser)

# find sensitivity given probabilities, nucleotide, and a filter threshold

# What is the accuracy of basecalling if I set a threshold a P=0.xx?
# Only applicable to unmodified IVT samples
find_accuracy <- function(calls_tsv, nt, threshold){
  
  #set up confusion matrix with all classes labels - true and false positive
  confusion_matrix = expand.grid(canonical_base = nt, class = c("TP", "FP"))
  
  calls_tsv %>% 
    filter(canonical_base == nt) %>% 
    #class nucleotides as PASS or FAIL by threshold
    mutate(prob_threshold = case_when(
      call_prob <= threshold ~ "FAIL",
      call_prob > threshold ~ "PASS"
    )) %>% 
    filter(prob_threshold == "PASS") %>% 
    group_by(canonical_base, call_code, prob_threshold) %>%
    reframe(code_total = n()) %>% 
    #class mod and canonical nucleotides into confusion matrix
    mutate(class = case_when(
      call_code == "-" & prob_threshold == "PASS" ~ "TP",
      call_code != "-" & prob_threshold == "PASS" ~ "FP"
    )) %>% 
    group_by(canonical_base, class) %>% 
    reframe(class_total = sum(code_total)) %>%
    #join with confusion matrix labels, if no rows exist then replace with zero
    full_join(confusion_matrix) %>% 
    replace_na(list(class_total = 0)) %>%
    group_by(canonical_base) %>% 
    #calc accuracy
    reframe(accuracy = class_total[class == "TP"] / sum(class_total) * 100) %>% 
    ungroup() %>% 
    pull(accuracy) %>% 
    round(digits = 5) %>% 
    suppressMessages()
}

# find threshold given probabilities, nucleotide and a target accuracy

# What threshold should I set if I want an accuracy of xx%?
find_optimum_threshold <- function(nt, calls_tsv, threshold_tsv, target_accuracy){
  
  # default threshold is the prob score of the 10th percentile
  default_threshold = threshold_tsv %>% 
    filter(base == nt & percentile == 10) %>% 
    pull(threshold) %>% 
    round(digits = 7) %>% as.numeric()
  
  default_accuracy = find_accuracy(calls_tsv, nt, default_threshold)
  
  # what is the accuracy when unfiltered
  raw_accuracy = find_accuracy(calls_tsv, nt, 0)
  
  # what is the maximum accuracy we can achieve with a maximum threshold? this is so we can say the target accuracy is too high and needs to be changed 
  max_threshold = 0.9999999
  max_accuracy = find_accuracy(calls_tsv, nt, max_threshold)
  
  
  # if target accuracy is higher than max accuracy we return the default and say try again
  if (max_accuracy < target_accuracy | default_threshold > max_threshold | is.na(default_accuracy)){
    final_threshold = max_threshold
    final_accuracy = max_accuracy
    
    # if default threshold is already good enough to achieve target accuracy AND the threshold is >0.7 (modkit recommended minimum), then we finish the search
  } else if (default_accuracy >= target_accuracy & default_threshold >= 0.7){
    final_threshold = default_threshold
    final_accuracy = default_accuracy
    
    # if default threshold is accurate enough but <0.7, then we set the threshold to be =0.7
  } else if (default_accuracy >= target_accuracy & default_threshold < 0.7){
    final_threshold = 0.7
    final_accuracy = find_accuracy(calls_tsv, nt, 0.7)
    
    # if default accuracy is lower than target, then we start the search
  } else if (default_threshold < max_threshold & default_accuracy < target_accuracy) {
    
    # set a range of threshold to calculate accuracy scores, starting with the default and the max 
    upperbound_threshold = max_threshold
    lowerbound_threshold = default_threshold
    
    # i is the counter fr the while loop but also the exponent for probability
    # starting with intervals of 0.001 (10^-3) between the lower and upper bound, we calculate accuracy for every 0.01 increment of threshold
    for (i in -3:-7){
        step = 10^i
        grid = seq(lowerbound_threshold, upperbound_threshold, step)
        accuracy_table = data.frame(prob = grid)
        accuracy_table = accuracy_table |> 
            mutate(accuracy = sapply(accuracy_table$prob, find_accuracy,
                                     calls_tsv = calls_tsv, nt = nt))

    # output all rows where accuracy is achieved, this keeps track of the running top candidates
    candidates = accuracy_table |> filter(accuracy >= target_accuracy)

    # if at least one target threshold is found, select lowest threshold
        if (nrow(candidates) > 0){
        best_candidate = candidates |> slice_min(prob)
        final_threshold = best_candidate |> pull(prob)
        final_accuracy = best_candidate |> pull(accuracy)
        
        # reset bounds to search in smaller steps
        lowerbound_threshold = max(default_threshold, final_threshold - step)
        upperbound_threshold = min(max_threshold, final_threshold + step)
        } else {
            # if there are no candidates found
            upperbound_threshold = max_threshold
            lowerbound_threshold = default_threshold
            }
    }
  }
  
  output_df <- data.frame(base = nt,
                          raw_accuracy = raw_accuracy, 
                          default_threshold = default_threshold, 
                          final_threshold = final_threshold, 
                          final_accuracy = final_accuracy)

   # # if final_threshold is the max at 0.9999999, we change it back to 1
   # mutate(final_threshold = case_match(final_threshold, 
   #                                     0.9999999 ~ 1,
   #                                     .default = final_threshold))
  
  return(output_df)
}

# run pipeline for ATCG, then output a single dataframe
run_for_single_chrom <- function(sample_chrom_dir){
  
  # for a given directory, find probs and threshold tsv
  calls_tsv = list.files(sample_chrom_dir, pattern = "extract_calls.tsv.gz", full.names = TRUE) %>% 
    read_table(show_col_types = FALSE,
      col_types = list(
      forward_read_position = "i",
      ref_position = "i",
      fw_soft_clipped_start = "i",
      fw_soft_clipped_end = "i",
      alignment_start = "i",
      alignment_end = "i",
      read_length = "i",
      base_qual = "i",
      flag = "i")) |> 
    as.data.frame() |> 
    filter(within_alignment == TRUE & ref_position != -1) |> 
    select(-c(forward_read_position, mod_strand, ref_mod_strand, fw_soft_clipped_start, fw_soft_clipped_end, alignment_start, alignment_end, read_length, base_qual, alignment_end, query_kmer, flag)) |> 
    group_by(chrom, call_code, canonical_base) |> 
    select(chrom, call_prob, call_code, canonical_base)

  threshold_tsv = list.files(sample_chrom_dir, pattern = "thresholds.tsv", full.names = TRUE) %>% 
    read_table(show_col_types = FALSE)

  threshold_tsv_nt = threshold_tsv |> distinct(base) |> pull()
  
  # find optimum threshold for all canonical nucleotides
  lapply(threshold_tsv_nt, find_optimum_threshold,
         calls_tsv = calls_tsv,
         threshold_tsv = threshold_tsv, 
         target_accuracy = opt$target_accuracy) %>% 
    # combine all nucleotide rows into one output
    bind_rows() %>%
    # add chrom column to the start
    add_column(chrom = basename(sample_chrom_dir), .before = 1) %>% 
    # combine all chrom AND ATCG rows into one output
    bind_rows()
}

run_pipeline_single_chrom <- function(sample_chrom_dir){

  flog.info(paste0("Analysing ", sample_chrom_dir))
  optimise_accuracy_df = run_for_single_chrom(sample_chrom_dir)
  
  if(is.na(opt$output_dir)){
    output_sample_dir = file.path(sample_chrom_dir)
  } else {
    # the user supplied an output dir
    output_sample_dir = file.path(opt$output_dir)
    #output_sample_dir = paste0(output_dir,"/",sample_name,"/")
    
    if (!dir.exists(output_sample_dir)) {
      dir.create(output_sample_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  write_tsv(optimise_accuracy_df, file = paste0(output_sample_dir,"/optimise_accuracy.tsv"))
  flog.info(paste0(sample_chrom_dir, " complete"))
}

#  DEPRECATED 
#  run pipeline for all subdirectories (chrom) of a single sample_dir, then output a single dataframe, and save output
run_pipeline <- function(sample_dir){
  
  sample_name = basename(sample_dir)
  flog.info(paste0("Analysing ", sample_name))
  
  chrom_dirs = list.dirs(sample_dir, recursive = FALSE)
  
  optimise_accuracy_df = lapply(chrom_dirs, run_for_single_chrom) %>% 
    # combine all chrom AND ATCG rows into one output
    bind_rows() %>% 
    # add sample name column to the start
    add_column(sample = sample_name, .before = 1)
  
  # find/create output dir with subdirectories
  if(is.na(opt$output_dir)){
    output_sample_dir = file.path(sample_dir)
  } else {
    # the user supplied an output dir
    output_sample_dir = file.path(opt$output_dir)
    #output_sample_dir = paste0(output_dir,"/",sample_name,"/")
    
    if (!dir.exists(output_sample_dir)) {
      dir.create(output_sample_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  write_tsv(optimise_accuracy_df, file = paste0(output_sample_dir,"/",sample_name,"_",opt$target_accuracy,".tsv"))
  flog.info(paste0(sample_name, " complete"))
  
}

# read input directory and run
if (file_test("-d", opt$input_dir) == TRUE){
  # it's a directory, recursively analyse all sample_dir
  
  # get a list of all sample_dir
  sample_chrom_dir = opt$input_dir
  #sample_dirs = list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
  
  flog.info("")
  flog.info("**** Starting sensitivity optimisation ****")
  
  mclapply(sample_chrom_dir, run_pipeline_single_chrom, mc.cores = opt$cores)
  
  #  # run pipeline for each sample
  #  final_output_df <- mclapply(sample_dir, run_pipeline, mc.cores = opt$cores) %>% 
  #    bind_rows()
  #  
  #  # save a big summary tsv for all samples
  #  if(is.na(opt$output_dir)){
  #    output_dir = file.path(input_dir)
  #  } else {
  #    # the user supplied an output dir
  #    output_dir = file.path(opt$output_dir)
  #  }
  #  
  #  write_tsv(final_output_df, file = paste0(output_dir,"probs_threshold_",opt$target_accuracy,".tsv"))
  
  flog.info("**** Analysis complete ****")
}
