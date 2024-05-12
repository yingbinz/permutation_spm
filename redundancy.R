
require(dplyr)

## get the sub and super associations

# judge if a pattern is a sub pattern, output is logical value, TRUE means is a subpattern
sub_matching <- function(remain_sub, remain_pattern, maxgap){
  if(length(remain_sub) > length(remain_pattern)){
    return(FALSE)
  } else if (length(remain_sub) == length(remain_pattern)){
    if (all(remain_sub == remain_pattern)){
      return(TRUE)
    } else{
      return(FALSE)
    }
  } else{
    if  (remain_sub[1] %in% remain_pattern[1:maxgap]){
      if(length(remain_sub) == 1){
        return(TRUE)
      } else {
        pos <- which(remain_pattern[1:maxgap] == remain_sub[1])[1]
        return(sub_matching(remain_sub = remain_sub[-1], remain_pattern = remain_pattern[-c(1:pos)], maxgap = maxgap))
      }
    } else{
      return(FALSE)
    } 
  }
}

########### get the pattern - sub pattern list,  #########
# get the sub and super associations
#   do while len > 2
#     for each of patterns with a len of longest 
#       judge which shorter patterns are its subpatterns 
#       record the subpatterns
#   longest -= 1
##################################################
get_pattern_sub_list <- function(sig_patterns, sig_p_lens, maxgap){
  # This is a function to identify the subpatterns of a pattern
  # input:
  #   sig_patterns: patterns that are significant
  #   sig_p_lens: the length of significant patterns, the order should be the same as sig_patterns 
  #   maxgap: the max gap
  # output:
  #   pattern_sub_list: a list of vectors, where the name of an vector is a pattern, and the elements are the pattern's subpatterns
  temp_longest <- sig_p_lens[which.max(sig_p_lens)]
  
  pattern_sub_list <- list()
  while (temp_longest > 2) {
    temp_patterns <- sig_patterns[which(sig_p_lens == temp_longest)]
    temp_subs <- sig_patterns[which(sig_p_lens == temp_longest-1)] #########sig_patterns[which(sig_p_lens < temp_longest-1)]
    for (pat in temp_patterns) {
      matched_subs <- c()
      pat_items <- unlist(strsplit(pat, ","))
      for (sub in temp_subs) {
        if(grepl(sub, pat)){
          matched_subs <- c(matched_subs, sub)
        } else{
          sub_items <- unlist(strsplit(sub, ","))
          matching <- sub_matching(remain_sub = sub_items, 
                                        remain_pattern = pat_items, 
                                        maxgap = maxgap)
          if(matching){
            matched_subs <- c(matched_subs, sub)
          } 
        }
      }
      pattern_sub_list[[pat]] <- matched_subs
    }
    temp_longest <- temp_longest - 1
  }
  
  return(pattern_sub_list)
}

# test the sig via permutation tests
compute_super_sub_diff <- function(raw_diff, permutated_df, pat, subs, redundancy){
  subs <- subs[! redundancy[subs]]
  diff_df <- permutated_df[,subs] - t(permutated_df[,pat])
  difference_p <- t(diff_df) > raw_diff[subs] 
  difference_p <- difference_p %>% t() %>% colMeans() %>% p.adjust(method = "BY") 
  difference_p <- difference_p <= 0.05
  names(difference_p) <- subs
  return(difference_p)
}

identify_redundancy <- function(sig_patterns, raw_mean_occurrences, permutation_df, maxgap){
  # This is a function to identify which patterns are redundant
  # input:
  #   sig_patterns: patterns that are significant
  #   raw_mean_occurrences: the mean occurrences of patterns in the raw data
  #   permutation_df: the permutated mean occurrences in each permutation, with column names as the pattern
  #   maxgap: the max gap
  # output:
  #   redundancy: a vector where the name of an item is a pattern, and the value is its redundancy
  #######################
  #   for each pattern that has subpatterns 
  #     test its difference with the subpatterns
  #       if the subpatterns are marked insignificant
  #         skip it
  #       else
  #         test its difference with its subpatterns#   
  #     marking insignificant subpatterns in redundancy
  ######################
  pat_redundant_sub_list <- c()
  sig_patterns <- sort(sig_patterns, decreasing = FALSE)
  sig_p_lens <- lapply(sig_patterns, function(x) length(strsplit(x, ",", fixed=TRUE, )[[1]])) %>% unlist()
  redundancy <- rep(FALSE, length(sig_patterns))
  names(redundancy) <- sig_patterns
  
  pattern_sub_list <- get_pattern_sub_list(sig_patterns=sig_patterns, sig_p_lens=sig_p_lens, maxgap=maxgap)
  
  
  for (pat in names(pattern_sub_list)) {
    subs <- pattern_sub_list[[pat]]
    raw_diff <- raw_mean_occurrences[subs] - raw_mean_occurrences[pat]
    permutated_sig <- compute_super_sub_diff(raw_diff = raw_diff,
                                          permutated_df = permutation_df, 
                                          pat = pat, 
                                          subs = subs,
                                          redundancy = redundancy)
    redundant_patterns <- names(permutated_sig)[!permutated_sig]
    if(length(redundant_patterns)){
      redundancy[redundant_patterns] <- TRUE
      pat_redundant_sub_list[[pat]] <- redundant_patterns
    }
  }
  
  
  
  return(list(redundancy, pat_redundant_sub_list))
}
