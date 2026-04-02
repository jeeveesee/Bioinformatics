# =========================
# Gibbs Sampler - tidyverse (final, tested structure)
# Names:
#   k -> kmer_len
#   t -> num_strings
#   profile -> profile_matrix
# =========================
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

# -------- helpers --------

# Accept list or character vector
.as_char_vec <- function(x) {
  if (is.list(x)) unlist(x, use.names = FALSE) else as.character(x)
}

# Convert motifs vector -> num_strings x kmer_len matrix (characters)
.motifs_mat <- function(motifs) {
  pieces   <- strsplit(motifs, "", fixed = TRUE)
  kmer_len <- length(pieces[[1]])
  matrix(unlist(pieces, use.names = FALSE),
         nrow = length(motifs),
         byrow = TRUE)
}

# Numerically stable log-sum-exp
.logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# -------- core pieces --------

# 1) Build profile_matrix with pseudocount = 1
tv_create_profile_matrix <- function(motifs) {
  nucs <- c("A","C","G","T")
  motifs <- .as_char_vec(motifs)
  num_strings <- length(motifs)
  M <- .motifs_mat(motifs)                 # num_strings x kmer_len
  kmer_len <- ncol(M)

  # counts per nucleotide per column -> 4 x kmer_len
  cnt <- vapply(nucs, function(n) colSums(M == n), numeric(kmer_len))
  cnt <- t(cnt)                             # rows A,C,G,T; cols 1..k

  profile_matrix <- (cnt + 1) / (num_strings + 4)   # pseudocounts
  rownames(profile_matrix) <- nucs
  profile_matrix
}

# 2) Weighted-random k-mer from one DNA string (uses paired indexing + log-space)
tv_most_probable_kmer <- function(dna, kmer_len, profile_matrix) {
  dna <- as.character(dna)[1]
  L <- nchar(dna)
  stops <- L - kmer_len + 1L
  if (stops < 1L) stop("kmer_len is longer than the sequence")

  # All candidate kmers
  kmers <- substring(dna, 1:stops, kmer_len:(kmer_len + stops - 1L))

  # Log-probability for each kmer via paired (row, col) picks
  logps <- vapply(kmers, function(kmer) {
    chars  <- strsplit(kmer, "", fixed = TRUE)[[1]]         # length kmer_len
    rix    <- match(chars, rownames(profile_matrix))        # A/C/G/T -> 1..4
    cix    <- seq_len(kmer_len)                             # 1..k
    vals   <- profile_matrix[cbind(rix, cix)]               # paired cells
    sum(log(vals))
  }, numeric(1))

  # Softmax to probs (stable)
  lse   <- .logsumexp(logps)
  probs <- exp(logps - lse)

  # Weighted sample
  sample(kmers, size = 1L, prob = probs)
}

# 3) Score motifs (lower is better)
tv_score_motifs <- function(motifs) {
  nucs <- c("A","C","G","T")
  motifs <- .as_char_vec(motifs)
  num_strings <- length(motifs)
  M <- .motifs_mat(motifs)
  kmer_len <- ncol(M)

  cnt <- vapply(nucs, function(n) colSums(M == n), numeric(kmer_len))
  cnt <- t(cnt)                                           # 4 x k
  max_per_col <- apply(cnt, 2, max)
  sum(num_strings - max_per_col)
}

# 4) Consensus from motifs
tv_consensus <- function(motifs) {
  nucs <- c("A","C","G","T")
  motifs <- .as_char_vec(motifs)
  M <- .motifs_mat(motifs)
  kmer_len <- ncol(M)

  cnt <- vapply(nucs, function(n) colSums(M == n), numeric(kmer_len))
  cnt <- t(cnt)                                           # 4 x k
  ixs <- apply(cnt, 2, which.max)
  paste0(nucs[ixs], collapse = "")
}

# 5) Main Gibbs sampler with multi-restarts
tv_gibbs_sampler <- function(dna_list, kmer_len, num_strings, N, num_restarts = 20L) {
  dna_list <- .as_char_vec(dna_list)
  Ls <- nchar(dna_list)

  best_overall_motifs <- NULL
  best_overall_score  <- Inf

  for (r in seq_len(num_restarts)) {
    # Random initialization
    starts <- map2_int(Ls, rep(kmer_len, length(Ls)),
                       \(L, kk) sample.int(L - kk + 1L, 1L))
    motifs <- map2_chr(dna_list, starts,
                       \(dna, s) str_sub(dna, s, s + kmer_len - 1L))

    best_motifs <- motifs
    best_score  <- tv_score_motifs(best_motifs)

    # Gibbs iterations
    for (iter in seq_len(N)) {
      i <- sample.int(num_strings, 1L)
      motifs_except_i <- motifs[-i]
      profile_matrix  <- tv_create_profile_matrix(motifs_except_i)
      motifs[i] <- tv_most_probable_kmer(dna_list[i], kmer_len, profile_matrix)

      current_score <- tv_score_motifs(motifs)
      if (current_score < best_score) {
        best_motifs <- motifs
        best_score  <- current_score
      }
    }

    if (best_score < best_overall_score) {
      best_overall_score  <- best_score
      best_overall_motifs <- best_motifs
    }
  }

  best_overall_motifs
}

# Convenience runner (same input layout that your Python example uses)
tv_run_gibbs_from_values <- function(dna_list, kmer_len, num_strings, N, num_restarts = 20L) {
  best <- tv_gibbs_sampler(dna_list, kmer_len, num_strings, N, num_restarts)
  cat("\nThe best motifs are: ", paste(best, collapse = " "), "\n", sep = "")
  cat("The consensus string is: ", tv_consensus(best), "\n", sep = "")
  invisible(list(best_motifs = best, consensus = tv_consensus(best)))
}

#############################################
#############################################
#                  TEST
#############################################
#############################################

dna_list <- list(
  'CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA',
  'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
  'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
  'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
  'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
)
kmer_len     <- 8
num_strings  <- 5
N            <- 100
num_restarts <- 20L

tv_run_gibbs_from_values(dna_list, kmer_len, num_strings, N, num_restarts)