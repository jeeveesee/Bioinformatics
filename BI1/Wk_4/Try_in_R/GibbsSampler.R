#################################################################################################
# TASK:
# Use Gibbs Sampler method to generate  profiles to
# find a set of motifs that are likely implanted by a common motif
# NOTE: YOU HAVE TO RUN THIS FROM THE BI1 FOLDER WITH THE COMMAND:
# python -m Wk_4.Wk4_2_GibbsSampler

# GibbsSampler(Dna, k, t, N)
#     randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
#     BestMotifs ← Motifs
#     for j ← 1 to N
#         i ← Random(t)
#         Profile ← profile matrix constructed from all strings in Motifs except for Motifi
#         Motifi ← Profile-randomly generated k-mer in the i-th sequence
#         if Score(Motifs) < Score(BestMotifs)
#             BestMotifs ← Motifs
#     return BestMotifs

#################################################################################################

library(glue)

# Main Gibbs Function
##################################################################################################

gibbs_sampler <- function(dna_list, k, t, N, num_restarts = 20) 
{
  best_overall_motifs <- NULL
  best_overall_score <- Inf

  for(restart in 1:num_restarts)
  {
    # Initialize random kmer motifs from each dna string
      motifs <- list()
      for(dna in dna_list)
      {
          # print(glue("Length of DNA is: {nchar(dna)}"))
          start_index <- sample(1:(nchar(dna) - k + 1), 1)
          print(start_index)
          motifs <- append(motifs, substr(dna, start_index, start_index + k -1))
          # print(unlist(motifs))
      }
      best_motifs <- motifs
    
    # Run Gibbs Sample with resubstitution N times
      for(n in 1:N)
      {
          i <- sample(1:t)
          motifs_except_i <- c(motifs[0:(i-1)], motifs[(i+1):(length(motifs))])
          profile <- create_profile_matrix(motifs_except_i)
          new_motif <- most_probable_kmer(dna_list[i], k, profile)
          motifs[i] <- new_motif

          if(score_motifs(motifs) < score_motifs(best_motifs))
          {
              best_motifs <- motifs  
          }
      }
      current_score <- score_motifs(best_motifs)
      if(current_score < best_overall_score)
      {
          best_overall_score <- current_score
          best_overall_motifs <- best_motifs
      }
    
    return(best_overall_motifs)
  }
}

# Helper Functions
##################################################################################################

  
# 1. Profile Matrix
###########################
  
  
gibbs_sampler(list('AAGTTCA', 'TTGGCCA', 'AGCCTCA'), 4, 3, 3, 1)
  