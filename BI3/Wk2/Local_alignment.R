# Local Alighment problem

# -------------------------------------------------------------------------
# 1. DEFINE THE PAM250 MATRIX
# -------------------------------------------------------------------------
create_pam250 <- function() {
  aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  
  # Full PAM250 raw data (20x20)
  values <- c(
     2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0, # A
    -2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2, # R
     0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3, -1,  1,  0, -4, -2, -2, # N
     0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2, # D
    -2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2, # C
     0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2, # Q
     0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2, # E
     1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1, # G
    -1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2, # H
    -1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4, # I
    -2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2, # L
    -1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2, # K
    -1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2, # M
    -3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1, # F
     1,  0, -1, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1, # P
     1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1, # S
     1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0, # T
    -6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6, # W
    -3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2, # Y
     0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4  # V
  )
  
  pam <- matrix(values, nrow = 20, ncol = 20, byrow = TRUE)
  rownames(pam) <- colnames(pam) <- aa
  return(pam)
}

# -------------------------------------------------------------------------
# 2. LOCAL ALIGNMENT FUNCTION (PURE R)
# -------------------------------------------------------------------------
solve_local_alignment <- function(v, w, sigma = 5) {
  pam250 <- create_pam250()
  
  v_seq <- strsplit(v, "")[[1]]
  w_seq <- strsplit(w, "")[[1]]
  n <- length(v_seq)
  m <- length(w_seq)
  
  # Scoring matrix S and Backtrack matrix
  S <- matrix(0, nrow = n + 1, ncol = m + 1)
  # 0: Stop, 1: Diag, 2: Up, 3: Left
  backtrack <- matrix(0, nrow = n + 1, ncol = m + 1)
  
  max_score <- 0
  max_pos <- c(1, 1)
  
  # Fill Matrix
  for (i in 2:(n + 1)) {
    v_char <- v_seq[i - 1]
    for (j in 2:(m + 1)) {
      w_char <- w_seq[j - 1]
      
      # Options: 0 (local restart), Match/Mismatch, Deletion, Insertion
      opts <- c(
        0, 
        S[i - 1, j - 1] + pam250[v_char, w_char], 
        S[i - 1, j] - sigma, 
        S[i, j - 1] - sigma
      )
      
      best_idx <- which.max(opts)

      cat(paste("i = ", i, "\n"))
      cat(paste("j = ", j, "\n"))
      cat(paste("v_char = ", v_char, "\n"))
      cat(paste("w_char = ", w_char, "\n"))
      cat(paste("best_idx = ", best_idx, "\n"))

      S[i, j] <- opts[best_idx]
      cat(paste("S[i,j] = ", S[i, j], "\n"))

      backtrack[i, j] <- best_idx - 1
      cat(paste("backtrack[i, j] = ", backtrack[i, j], "\n"))

      if (S[i, j] > max_score) {
        max_score <- S[i, j]
        max_pos <- c(i, j)
      }
    }
  }
  
  # Traceback
  res_v <- character(n + m)
  res_w <- character(n + m)
  curr_i <- max_pos[1]
  curr_j <- max_pos[2]
  k <- n + m
  
  while (curr_i > 1 || curr_j > 1) {
    if (S[curr_i, curr_j] == 0) break
    
    dir <- backtrack[curr_i, curr_j]
    if (dir == 1) { # Diagonal
      res_v[k] <- v_seq[curr_i - 1]
      res_w[k] <- w_seq[curr_j - 1]
      curr_i <- curr_i - 1; curr_j <- curr_j - 1
    } else if (dir == 2) { # Up
      res_v[k] <- v_seq[curr_i - 1]
      res_w[k] <- "-"
      curr_i <- curr_i - 1
    } else if (dir == 3) { # Left
      res_v[k] <- "-"
      res_w[k] <- w_seq[curr_j - 1]
      curr_j <- curr_j - 1
    }
    k <- k - 1
  }
  
  # Final formatting
  cat(max_score, "\n")
  cat(paste(res_v[(k+1):(n+m)], collapse = ""), "\n")
  cat(paste(res_w[(k+1):(n+m)], collapse = ""), "\n")
}

# -------------------------------------------------------------------------
# 3. TEST WITH SAMPLE DATA
# -------------------------------------------------------------------------
v_input <- "MEANLY"
w_input <- "PENALTY"

solve_local_alignment(v_input, w_input)