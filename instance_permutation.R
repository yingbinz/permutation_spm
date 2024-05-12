require(dplyr)
require(purrr)
require(arulesSequences)
# not 100% accurate when gap>2, but error rate lower than 0.013%, and bias is no more than 1 occurrences even a sequence contains hundreds of events.
# issues in compute_EID_in_shared and find_overlapping

### get frequent sequential patterns ####
getFreqSeq <- function(data, support, maxgap, ID_col, EID_col, event_col, remove_single =TRUE,dataname=""){
  # if(nrow(data)>99999){
  #   data[,EID_col] <- format(data[,EID_col], scientific = FALSE)
  # }
  temp_file_name <- paste0("formated_event_data",dataname,".txt")
  data %>% dplyr::select(!!sym(ID_col), !!sym(EID_col), !!sym(event_col)) %>% 
    arrange(!!sym(ID_col), !!sym(EID_col)) %>% 
    write.table(., temp_file_name, sep=";", 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  event_data_matrix <- read_baskets(temp_file_name, 
                                    sep = ";", 
                                    info =c("sequenceID","eventID"))
  
  freq_seq <- cspade(event_data_matrix, parameter = list(support = support, 
                                                         maxsize = 1, maxgap = maxgap), 
                     control = list(verbose = F))
  # freq_seq <- ruleInduction(freq_seq, confidence = 0.001,
  #                           control = list(verbose = FALSE))
  freq_seq <- as(freq_seq, "data.frame") %>% arrange(desc(support)) 
  colnames(freq_seq)[1] <- "sequence"
  freq_seq$length <- map_dbl(freq_seq$sequence,
                             function(x) strsplit(as.character(x), ',', fixed = TRUE) %>%
                               unlist() %>% length())
  
  if(remove_single){freq_seq <- subset(freq_seq, length>1)}
  
  return(freq_seq)
}

### preprocess data ####
preprocess <- function(data, ID_col, event_col, EID_col, gap){
  for (gap_i in 1:gap) {
    next_ID_n <- paste('next',gap_i,ID_col,sep='_')
    next_event_n <- paste('next',gap_i,'event',sep='_')
    next_eventID_n <- paste('next',gap_i,'EID',sep='_')
    data[, next_ID_n] <- c(unlist(data[-c(1:gap_i), ID_col]), rep("No ID",gap_i))
    data[,next_event_n] <- c(unlist(data[-c(1:gap_i), event_col]), rep(-9, gap_i))
    data[, next_eventID_n] <- c(unlist(data[-c(1:gap_i), EID_col]), rep(-9,gap_i))
    
    unsuccessive <- (data[,EID_col] != data[,next_eventID_n] - gap_i) |
      (data[,ID_col] != data[,next_ID_n])
    
    data[unsuccessive, next_event_n] <-  -9
    data[unsuccessive, next_eventID_n] <-  -9
  }
  
  return(data)
}

### functions for counting pattern occurrences ####

searching_used_ids <- function(used_event_id_c, re2_matched, 
                               event_n, temp_used_event_id_c,
                               num_duplicated_EID_c,
                               re2_c, pattern_len_m1, used_event_id, 
                               next_eventID_c, match_eventID_c, EID_col){
  
  found_nonoverlapping <- FALSE
  can_searching_deeper <- (event_n<pattern_len_m1) && any(num_duplicated_EID_c[-c(1:event_n)]>0)
  
  candidate_event_id <- re2_matched[next_eventID_c][
    as.logical(re2_matched[match_eventID_c])]
  
  for (used_event_id in candidate_event_id) {
    if(!(used_event_id %in% used_event_id_c[[event_n]])){
      temp_used_event_id_c <- c(temp_used_event_id_c, used_event_id)
      if(can_searching_deeper){
        
        next_re2_matched <- re2_c[[event_n+1]] 
        next_re2_matched <- next_re2_matched[next_re2_matched[,EID_col]==used_event_id, ]
        searching_result <- searching_used_ids(used_event_id_c, next_re2_matched, 
                                               event_n+1, temp_used_event_id_c,
                                               num_duplicated_EID_c,
                                               re2_c, pattern_len_m1, used_event_id,
                                               next_eventID_c, match_eventID_c, EID_col)
        
        found_nonoverlapping <- searching_result[[2]]
        
        if(found_nonoverlapping){ 
          temp_used_event_id_c <- searching_result[[1]]
          break 
        }
        
      } else{
        found_nonoverlapping <- TRUE
        break
      }
      ## return which event has been used
      # if( (used_event_id == candidate_event_id[length(candidate_event_id)]) &
      #     ){
      #   next
      # } else{
      temp_used_event_id_c <- temp_used_event_id_c[-length(temp_used_event_id_c)]
      # }
      
    }
    
  }
  
  return(list(temp_used_event_id_c, found_nonoverlapping))  
}

patterns_to_events <- function(patterns){

  events_in_patterns <- patterns %>% 
    strsplit(",", fixed=TRUE) %>%
    lapply(., as.numeric) 
  
  patterns <- map_chr(events_in_patterns, function(x) paste(event_dic[x],collapse = ","))
  
  return(list(patterns=patterns, events_in_patterns=events_in_patterns))
}


find_unique <- function(temp_re2_c,
                        EID_col,
                        next_eventID_c, match_eventID_c){
  final_EID_df <- cbind(next_eventID=c(temp_re2_c[,next_eventID_c]), matching=c(temp_re2_c[,match_eventID_c]))
  
  final_EID_df <- final_EID_df[final_EID_df[,"matching"]==1,]
  if(!is.matrix(final_EID_df)){
    final_EID_df <- t(as.matrix(final_EID_df))
  }   
  unique_EIDs <- unique(final_EID_df[,"next_eventID"])
  return(unique_EIDs)
}

find_overlapping <- function(temp_re2_c,
                             EID_col,
                             next_eventID_c, match_eventID_c,gap){
  # overlapping has been checked before. if a eid occur in a row where only one matching, it is impossible a possible overlapping
  # so only checked duplicate eid only in rows with #matching > 1, or only in rows with #matching == 1
  more_matching <- temp_re2_c[,match_eventID_c]
  if(is.matrix(more_matching)){
    more_matching <- rowSums(more_matching) > 1
  }  else{
    more_matching <- sum(more_matching)> 1
  }
  
  final_EID_df <- cbind(next_eventID=c(temp_re2_c[more_matching,next_eventID_c]), 
                        matching=c(temp_re2_c[more_matching,match_eventID_c]))
  final_EID_df <- final_EID_df[final_EID_df[,"matching"]==1,]
  if(!is.matrix(final_EID_df)){
    final_EID_df <- t(as.matrix(final_EID_df))
  }    
  duplicate_rows <- duplicated(final_EID_df[,"next_eventID"])
  possible_overlapping_1 <- unique(final_EID_df[duplicate_rows,"next_eventID"])
  
  final_EID_df <- cbind(next_eventID=c(temp_re2_c[!more_matching,next_eventID_c]), 
                        matching=c(temp_re2_c[!more_matching,match_eventID_c]))
  final_EID_df <- final_EID_df[final_EID_df[,"matching"]==1,]
  if(!is.matrix(final_EID_df)){
    final_EID_df <- t(as.matrix(final_EID_df))
  }    
  duplicate_rows <- duplicated(final_EID_df[,"next_eventID"])
  possible_overlapping_2 <- unique(final_EID_df[duplicate_rows,"next_eventID"])
  
  return(unique(c(possible_overlapping_2, possible_overlapping_1)))
}

find_ids_with_possible_overlapping <- function(re2_c_plen_2, next_eventID_c,
                                               possible_overlapping_EIDs, gap, EID_col){
  # two possible overlapping. First, the next event is in next_possible_overlapping_EIDs
  temp_overlapping <- matrix(re2_c_plen_2[,next_eventID_c] %in% possible_overlapping_EIDs, ncol=gap, byrow=FALSE)
  temp_overlapping <- apply(temp_overlapping, 1, any)
  temp_overlapping <- re2_c_plen_2[temp_overlapping, EID_col]
  return(temp_overlapping)
}

compute_EID_in_shared <- function(re2_c_plen_2, match_eventID_c, 
                                  next_eventID_c, shared_EIDs, gap){
  
  in_shared <- matrix(re2_c_plen_2[,next_eventID_c] %in% shared_EIDs, ncol=gap, byrow=FALSE)
  re2_c_plen_2[,match_eventID_c] <- in_shared
  in_shared <- c(apply(in_shared, 1, any))
  ## some EIDs could be used by multiple rows, let rows with only one matching use these EIDS firstly. Thus create 
  ## shared_EIDs_partial that are not used by rows with only one matching. Rows with more than one matching
  ## need to have EIDs in shared_EIDs_partial to be counted as in shared EIDs
  shared_EIDs_partial <- shared_EIDs
  for (i_gap in 1:(gap-1)) {
    more_matching <- re2_c_plen_2[,match_eventID_c]
    if(is.matrix(more_matching)){
      more_matching <- rowSums(more_matching) > 1
    }  else{
      more_matching <- sum(more_matching)> 1
    }
    shared_EIDs_partial <- shared_EIDs_partial[! shared_EIDs_partial %in% c(re2_c_plen_2[!more_matching,next_eventID_c])]
    in_shared_partial <- matrix(re2_c_plen_2[more_matching,next_eventID_c] %in% shared_EIDs_partial, ncol=gap, byrow=FALSE)
    re2_c_plen_2[more_matching, match_eventID_c] <- in_shared_partial
    
    in_shared_partial <- c(apply(in_shared_partial, 1, any))
    
    in_shared[more_matching] <- in_shared_partial
  }
  
  re2_c_plen_2 <- re2_c_plen_2[in_shared,]
  if(!is.matrix(re2_c_plen_2)){
    re2_c_plen_2 <- t(as.matrix(re2_c_plen_2))
  }
  return(re2_c_plen_2)
}
event_matching <- function(event_position, events_in_pattern, 
                           re2_c_list=NULL,unchecked_re2_c_list=NULL,
                           re0=NULL, 
                           rawdata, gap,
                           event_col, EID_col, pattern_len,
                           next_eventID_c, match_eventID_c, pattern_len_m1){
  ## event position = 1####
  if (event_position == 1){
    ### pattern_len = 2####
    if(pattern_len == 2){
      #from rawdata, extract all rows matching the first event of pattern: re2
      re2 <- subset(rawdata, rawdata[,event_col]== events_in_pattern[event_position])
      # whether the other events of pattern occur
      # event_position=event_position+1
      # re0=re2
      
      re2_matching <- event_matching(event_position=event_position+1, 
                                     events_in_pattern, 
                                     re2_c_list=NULL,unchecked_re2_c_list=NULL,
                                     re0=re2, 
                                     rawdata, gap,
                                     event_col, EID_col, pattern_len,
                                     next_eventID_c, match_eventID_c, pattern_len_m1)
      re2 <- cbind(re2, matching =re2_matching[[1]])
      
      re2_c <- re2_matching[[2]]
      duplicated_EID_c <- re2_matching[[3]]
      possible_overlapping_EIDs <- re2_matching[[4]]
      
    } else{
      
      ## create re2_c and duplicated_EID_c ####
      re2_c_and_post <- c(re2_c_list[paste0(events_in_pattern[1:(pattern_len-1)], collapse =",")],
                          unchecked_re2_c_list[paste0(events_in_pattern[(pattern_len-1):pattern_len], collapse =",")])
      
      duplicated_EID_c <- list()
      re2_c <- re2_c_and_post[[1]][[1]]
      unique_EIDs_list <- re2_c_and_post[[1]][[2]]
      re2_c_plen_2 <-  re2_c_and_post[[2]][[1]][[1]]
      shared_EIDs <- intersect(unique_EIDs_list[[pattern_len-2]], re2_c_plen_2[,EID_col])
      ### pattern_len 3 ####
      #### update the last item ####
      re2_c_plen_2 <- re2_c_plen_2[re2_c_plen_2[,EID_col] %in% shared_EIDs, ]
      if(!is.matrix(re2_c_plen_2)){
        re2_c_plen_2 <- t(as.matrix(re2_c_plen_2))
      }
      re2_c[[pattern_len-1]] <- re2_c_plen_2
      ##### cal the possible_overlapping_EIDs ####
      if(gap>1){
        possible_overlapping_EIDs <- cbind(next_eventID=c(re2_c_plen_2[,next_eventID_c]), 
                                           matching=c(re2_c_plen_2[,match_eventID_c]))
        
        possible_overlapping_EIDs <- possible_overlapping_EIDs[possible_overlapping_EIDs[,"matching"]==1,]
        if(!is.matrix(possible_overlapping_EIDs)){
          possible_overlapping_EIDs <- t(as.matrix(possible_overlapping_EIDs))
        }        
        duplicate_rows <- duplicated(possible_overlapping_EIDs[,"next_eventID"])
        possible_overlapping_EIDs <- possible_overlapping_EIDs[duplicate_rows,"next_eventID"]
        duplicated_EID_c[[pattern_len-1]] <- possible_overlapping_EIDs
        
        # two possible overlapping. First, the next event is in next_possible_overlapping_EIDs
        temp_overlapping <- find_ids_with_possible_overlapping(re2_c_plen_2, next_eventID_c,
                                                               possible_overlapping_EIDs, gap, EID_col)
      }
      
      #### update the last second item ####
      re2_c_plen_2 <- re2_c[[pattern_len-2]]
      
      if(gap>1){
        re2_c_plen_2 <- compute_EID_in_shared(re2_c_plen_2, match_eventID_c, next_eventID_c, shared_EIDs, gap)
      } else{
        in_shared <- re2_c_plen_2[,next_eventID_c] %in% shared_EIDs
        re2_c_plen_2 <- re2_c_plen_2[in_shared,]
        if(!is.matrix(re2_c_plen_2)){
          re2_c_plen_2 <- t(as.matrix(re2_c_plen_2))
        }
      }
      
      re2_c[[pattern_len-2]] <- re2_c_plen_2
      ##### cal the possible_overlapping_EIDs ####
      if(gap>1){
        
        # Second, the current event is found possible overlapping by find_overlapping        
        possible_overlapping_EIDs <- find_overlapping(temp_re2_c = re2_c_plen_2,
                                                      EID_col = EID_col, next_eventID_c=next_eventID_c, 
                                                      match_eventID_c = match_eventID_c)
        possible_overlapping_EIDs <- possible_overlapping_EIDs[possible_overlapping_EIDs %in% shared_EIDs]
        
        # combine two types of overlapping
        possible_overlapping_EIDs <- unique(c(possible_overlapping_EIDs, temp_overlapping))
        duplicated_EID_c[[pattern_len-2]] <- possible_overlapping_EIDs
        
        # two possible overlapping. First, the next event is in next_possible_overlapping_EIDs
        temp_overlapping <- find_ids_with_possible_overlapping(re2_c_plen_2, next_eventID_c,
                                                               possible_overlapping_EIDs, gap, EID_col)
      } 
      
      ### pattern_len > 3 ####
      if(pattern_len>3){
        for (event_n in (pattern_len-3):1) {
          shared_EIDs <- intersect(unique_EIDs_list[[event_n]], re2_c_plen_2[,EID_col])
          re2_c_plen_2 <- re2_c[[event_n]]
          if(gap>1){
            re2_c_plen_2 <- compute_EID_in_shared(re2_c_plen_2, match_eventID_c, next_eventID_c, shared_EIDs, gap)
          } else{
            in_shared <- re2_c_plen_2[,next_eventID_c] %in% shared_EIDs
            re2_c_plen_2 <- re2_c_plen_2[in_shared,]
            if(!is.matrix(re2_c_plen_2)){
              re2_c_plen_2 <- t(as.matrix(re2_c_plen_2))
            }
          }
          ### update the last second item ####
          re2_c[[event_n]] <- re2_c_plen_2
          #### cal the possible_overlapping_EIDs ####
          if(gap>1){
            # Second, the current event is found possible overlapping by find_overlapping
            possible_overlapping_EIDs <- find_overlapping(temp_re2_c = re2_c_plen_2,
                                                          EID_col = EID_col, next_eventID_c=next_eventID_c, 
                                                          match_eventID_c = match_eventID_c)
            possible_overlapping_EIDs <- possible_overlapping_EIDs[possible_overlapping_EIDs %in% shared_EIDs]
            
            # combine two types of overlapping
            possible_overlapping_EIDs <- unique(c(possible_overlapping_EIDs, temp_overlapping))
            
            duplicated_EID_c[[event_n]] <- possible_overlapping_EIDs
            
            # two possible overlapping. First, the next event is in next_possible_overlapping_EIDs
            temp_overlapping <- find_ids_with_possible_overlapping(re2_c_plen_2, next_eventID_c,
                                                                   possible_overlapping_EIDs, gap, EID_col)
          }
        }
      }
      if(gap>1){
        if(length(duplicated_EID_c[[1]])==0){
          duplicated_EID_c <- list()
        }
      } 
      re2 <- re2_c[[1]]
      re2 <- re2[,-which(colnames(re2) %in% match_eventID_c)]
      if(!is.matrix(re2)){
        re2 <- t(as.matrix(re2))
      }
    }  
    ###### check overlapping########        
    if((gap > 1)& (length(unlist(duplicated_EID_c)) > 0) ){
      num_duplicated_EID_c <- sapply(duplicated_EID_c, length)
      
      re2_matched <- re2_c[[1]]
      re2_matched <- cbind(re2_matched, possible_overlapping = 0)
      
      for (gap_i in 1:gap) {
        next_eventID_n <- paste('next',gap_i,'EID',sep='_')
        re2_matched[re2_matched[, next_eventID_n] %in% possible_overlapping_EIDs, 
                    'possible_overlapping'] <- 1
      }
      
      re2_matched <- subset(re2_matched, re2_matched[, 'possible_overlapping'] & 
                              re2_matched[, 'matching'])
      used_event_id_c <- rep(list(NULL), pattern_len_m1)
      found_c <- c()
      
      for (row_i in 1:nrow(re2_matched)) {
        temp_used_event_id_c <- c()
        found <- FALSE
        
        searching_result <- searching_used_ids(used_event_id_c, re2_matched[row_i,],
                                               event_n=1, temp_used_event_id_c,
                                               num_duplicated_EID_c,
                                               re2_c, pattern_len_m1, used_event_id=NULL,
                                               next_eventID_c, match_eventID_c, EID_col)
        
        temp_used_event_id_c <- searching_result[[1]]
        found <- searching_result[[2]]
        
        # if found non overlapping occurrences, update used_event_id_c
        found_c <- c(found_c, found)
        if(found){
          for (event_n in 1:pattern_len_m1) {
            temp_used_event_id <- temp_used_event_id_c[event_n]
            used_event_id_c[[event_n]] <- c(used_event_id_c[[event_n]], temp_used_event_id)
          }
        }
      }
      
      overlapping_id_c <- re2_matched[!found_c, EID_col]
      if(pattern_len==2){
        unchecked_re2 <- re2
      }
      
      re2[re2[,EID_col] %in% overlapping_id_c, 'matching'] <- 0
      
    }
    
  }  else{
    re2 <- re0
    re2 <- cbind(re2, matching=0)
    # re2[,'matching'] <- FALSE
    EID_c <- c()
    for (gap_i in 1:gap) {
      match_eventID_n <- paste('match', gap_i, sep="_")
      next_event_n <- paste('next',gap_i,'event',sep='_')
      next_eventID_n <- paste('next',gap_i,'EID',sep='_')
      # whether at least one of the next n rows has the event_position event
      temp_matching <- re2[,next_event_n] == events_in_pattern[event_position]
      # re2 <- re2 %>% mutate(matching = if_else(!!sym(match_eventID_n), TRUE, matching))
      re2[ temp_matching, 'matching'] <- 1
      
      EID_c <- c(EID_c, re2[temp_matching, next_eventID_n])
      
      # cbind at the end of a iter because re2 is matrix, cbind will cause temp_matching to 0/1
      re2 <- cbind(re2, temp_matching)
      colnames(re2)[ncol(re2)] <- match_eventID_n
    }
    duplicated_EID <- EID_c[duplicated(EID_c)]
    possible_overlapping_EID_c <- duplicated_EID
    duplicated_EID_c <- list(duplicated_EID)
    re2_c <- list(re2)
  }
  
  matching_EID_c <- re2[as.logical(re2[,'matching']),EID_col]
  
  #for each row in re0, whether at least one of the next n EID_col in matching_EID_c
  re0_matching <- re0[,EID_col] %in% matching_EID_c
  
  if(event_position == 1){
    
    re2_c_new <- re2_c
    temp_re2_c <- re2_c_new[[1]][as.logical(re2[,'matching']),]
    
    if(!is.matrix(temp_re2_c)){
      temp_re2_c <- t(as.matrix(temp_re2_c))
    }
    
    re2_c_new[[1]] <- temp_re2_c
    if(gap>1){    
      unique_EIDs <- find_unique(temp_re2_c = temp_re2_c,
                                 EID_col = EID_col, next_eventID_c=next_eventID_c, 
                                 match_eventID_c = match_eventID_c )
    } else{
      unique_EIDs <- temp_re2_c[,next_eventID_c]
    }
    
    unique_EIDs_list <- list(unique_EIDs)
    
    if(pattern_len==2){
      if(gap>1){
        if(!exists("unchecked_re2")){
          unchecked_re2 <- re2
        }
        unchecked_re2_c_new <- re2_c
        unchecked_re2_c <- unchecked_re2_c_new[[1]][as.logical(unchecked_re2[,'matching']),]
        unchecked_re2_c_new[[1]] <- unchecked_re2_c
        
        unchecked_unique_EIDs <- find_unique(temp_re2_c = unchecked_re2_c,
                                             EID_col = EID_col, next_eventID_c=next_eventID_c, 
                                             match_eventID_c = match_eventID_c )
        unchecked_unique_EIDs_list <- list(unchecked_unique_EIDs)        
      } else{
        unchecked_re2_c_new <- re2_c_new
        unchecked_unique_EIDs_list <- unique_EIDs_list
      }
    }
    
    if(pattern_len>2){
      if(gap>1){
        for (event_n in 2:pattern_len_m1) {
          temp_re2_c <- re2_c_new[[event_n]]
          if(!is.matrix(temp_re2_c)){
            temp_re2_c <- t(as.matrix(temp_re2_c))
          }
          temp_re2_c <- temp_re2_c[temp_re2_c[,EID_col] %in% unique_EIDs, ]
          if(!is.matrix(temp_re2_c)){
            temp_re2_c <- t(as.matrix(temp_re2_c))
          }
          re2_c_new[[event_n]] <- temp_re2_c
          
          unique_EIDs <- find_unique(temp_re2_c = temp_re2_c,
                                     EID_col = EID_col, next_eventID_c=next_eventID_c, 
                                     match_eventID_c = match_eventID_c )
          unique_EIDs_list[[event_n]] <- unique_EIDs
        }
      } else{
        for (event_n in 2:pattern_len_m1) {
          temp_re2_c <- re2_c_new[[event_n]]
          if(!is.matrix(temp_re2_c)){
            temp_re2_c <- t(as.matrix(temp_re2_c))
          }
          unique_EIDs <- temp_re2_c[,next_eventID_c]
          unique_EIDs_list[[event_n]] <- unique_EIDs
        }        
      }
      return (list(re0_matching, list(re2_c_new, unique_EIDs_list))) 
      
    } else{
      return (list(re0_matching, list(re2_c_new, unique_EIDs_list), list(unchecked_re2_c_new, unchecked_unique_EIDs_list))) 
    }
    
  } else{
    return (list(re0_matching, re2_c, duplicated_EID_c, possible_overlapping_EID_c)) 
  }
}


count_instances <- function(events_in_patterns, event_dic, data, ID_col, event_col, EID_col, gap){
  data <- data %>% preprocess(., ID_col = ID_col, event_col = event_col,
                              EID_col = EID_col, gap = gap)
  
  
  data_ids <- data[,ID_col]
  data <- dplyr::select(data, -contains(ID_col)) %>% as.matrix()
  
  patterns_lens <- sapply(events_in_patterns, length)
  
  next_eventID_c <- paste('next',1:gap,'EID',sep='_')
  match_eventID_c <- paste('match', 1:gap, sep="_")
  ## length 2 patterns####  
  temp_data <- lapply(events_in_patterns[which(patterns_lens==2)], 
                      function(ep) event_matching(event_position=1, events_in_pattern=ep,
                                                  re2_c_list=NULL,unchecked_re2_c_list=NULL,
                                                  re0=data, rawdata=data, gap=gap,
                                                  event_col, EID_col, pattern_len = 2,
                                                  next_eventID_c, match_eventID_c, pattern_len_m1=1))
  sorted_events_in_patterns <- events_in_patterns[which(patterns_lens==2)]
  temp_names <- sapply(events_in_patterns[which(patterns_lens==2)], function(x) paste0(x, collapse =","))
  instance_vectors <- lapply(temp_data, function(x) x[[1]])
  names(instance_vectors) <- temp_names
  instance_vectors <- bind_cols(instance_vectors)
  instance_vectors[,ID_col] <- data_ids
  instances <- instance_vectors %>%
    group_by(!!sym(ID_col)) %>%
    summarise_all(sum, na.rm = TRUE) %>% ungroup()
  
  re2_c_list <- lapply(temp_data, function(x) x[[2]])
  names(re2_c_list) <- temp_names
  
  unchecked_re2_c_list <- lapply(temp_data, function(x) x[[3]])
  names(unchecked_re2_c_list) <- temp_names
  
  max_len <- max(patterns_lens)
  ## length 3 patterns####  
  if(max_len > 2){
    for (pattern_len in 3:max_len) {
      remove(temp_data)
      
      temp_eps <- events_in_patterns[which(patterns_lens==pattern_len)]
      
      temp_data <- lapply(temp_eps, 
                          function(ep) event_matching(1, events_in_pattern=ep,
                                                      re2_c_list=re2_c_list,unchecked_re2_c_list=unchecked_re2_c_list,
                                                      re0=data, rawdata=data, gap=gap,
                                                      event_col, EID_col, pattern_len = pattern_len,
                                                      next_eventID_c, match_eventID_c, pattern_len_m1=pattern_len-1))
      
      sorted_events_in_patterns <- c(sorted_events_in_patterns, temp_eps)
      temp_names <- sapply(temp_eps, function(x) paste0(x, collapse =","))
      
      instance_vectors <- lapply(temp_data, function(x) x[[1]])
      names(instance_vectors) <- temp_names
      instance_vectors <- bind_cols(instance_vectors)
      instance_vectors[,ID_col] <- data_ids
      instance_vectors <- instance_vectors %>%
        group_by(!!sym(ID_col)) %>%
        summarise_all(sum, na.rm = TRUE) %>% ungroup()
      instances <- left_join(instances, instance_vectors, by = ID_col)
      temp_re2_c_list <- lapply(temp_data, function(x) x[[2]])
      names(temp_re2_c_list) <- temp_names
      re2_c_list <- c(re2_c_list, temp_re2_c_list)
    }    
  }
  
  patterns <- map_chr(sorted_events_in_patterns, function(x) paste(event_dic[x],collapse = ","))
  
  colnames(instances) <- c(ID_col, patterns)
  return(instances)
}



permutation <- function(event_seq_base, ID_col, event_col, seed=NULL, weighted=FALSE, w = NULL){
  set.seed(seed)
  permutated_data  <- event_seq_base %>%
    group_by_at({{ID_col}}) %>%
    mutate(!!event_col :=sample(!!sym(event_col))) %>%
    ungroup()

  return(permutated_data)
}

permutating_and_counting <- function(data, gap , ID_col , event_col , EID_col, 
                                     events_in_patterns , event_dic ,  seed=NULL){
  set.seed(seed = seed)
  permutated_data  <- data %>% permutation(data, ID_col = ID_col, event_col = event_col,seed=seed)
  
  temp_instances <- count_instances(events_in_patterns = events_in_patterns,
                                    event_dic = event_dic,
                                    data = permutated_data,
                                    ID_col = ID_col, event_col = event_col, EID_col = EID_col,
                                    gap = maxgap)
  temp_instances <- colMeans(as.matrix(temp_instances[,-1]))
  return(temp_instances)
}