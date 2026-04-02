# MULTIPLE PLACEHOLDERS:
# MAGRITTR: dat_with_super_long_name %>% {plot(.$x, .$y, cex = 3, lwd = 5)}
# BASE R EQUIVALENT: dat_with_super_long_name |> 
# (\(.) plot(.$x, .$y, cex = 3, lwd = 5))()

# Branch-and-Bound Cyclopeptide Sequencing Algorithm in R

# Amino acid mass table (integer masses)
AMINO_ACID_MASS <- c(
  G = 57, A = 71, S = 87, P = 97, V = 99, T = 101, C = 103, I = 113, L = 113,
  N = 114, D = 115, K = 128, Q = 128, E = 129, M = 131, H = 137, F = 147,
  R = 156, Y = 163, W = 186
)

# Calculate the mass of a peptide (vector of amino acid masses)
Mass <- function(peptide) 
        {
          if (length(peptide) == 0) return(0)
          return(sum(peptide))
        }

# Get parent mass (largest value in spectrum)
ParentMass <- function(spectrum) 
{
  return(max(spectrum))
}

# Generate theoretical linear spectrum
LinearSpectrum <- function(peptide) 
{
  if (length(peptide) == 0) return(c(0))
  
  prefix_mass <- c(0, cumsum(peptide))
  linear_spectrum <- c(0)
  
  for (i in 1:length(peptide)) 
    {
    for (j in (i+1):(length(peptide)+1)) 
      {
      linear_spectrum <- c(linear_spectrum, prefix_mass[j] - prefix_mass[i])
      }
    }
  
  return(sort(linear_spectrum))
}

# Generate theoretical cyclic spectrum
Cyclospectrum <- function(peptide) 
{
  if (length(peptide) == 0) return(c(0))
  
  prefix_mass <- c(0, cumsum(peptide))
  peptide_mass <- prefix_mass[length(prefix_mass)]
  cyclic_spectrum <- c(0, peptide_mass)
  
  for (i in 1:length(peptide)) 
    {
    for (j in (i+1):(length(peptide)+1)) 
      {
      mass <- prefix_mass[j] - prefix_mass[i]
      cyclic_spectrum <- c(cyclic_spectrum, mass)
      
      if (i > 1 && j <= length(peptide)) 
        {
        cyclic_spectrum <- c(cyclic_spectrum, peptide_mass - mass)
        }
      }
    }
  
  return(sort(cyclic_spectrum))
}

# Check if peptide is consistent with spectrum
IsConsistent <- function(peptide, spectrum) 
{
  linear_spec <- LinearSpectrum(peptide)
  
  # Create frequency tables
  spectrum_freq <- table(spectrum)
  linear_freq <- table(linear_spec)
  
  # Check if all masses in linear spectrum are in experimental spectrum
  # with sufficient frequency
  for (mass in names(linear_freq)) 
    {
    if (!(mass %in% names(spectrum_freq))) 
    {
      return(FALSE)
    }
    if (linear_freq[mass] > spectrum_freq[mass]) 
    {
      return(FALSE)
    }
  }
  
  return(TRUE)
}

# Expand candidate peptides by adding each possible amino acid mass
Expand <- function(peptides, masses = AMINO_ACID_MASS) 
{
  if (length(peptides) == 0) 
  {
    return(lapply(unique(masses), function(m) c(m)))
  }
  
  expanded <- list()
  unique_masses <- unique(as.vector(masses))
  
  for (peptide in peptides) 
  {
    for (mass in unique_masses) 
    {
      expanded[[length(expanded) + 1]] <- c(peptide, mass)
    }
  }
  
  return(expanded)
}

# Compare two spectra for equality
SpectraEqual <- function(spec1, spec2) 
{
  if (length(spec1) != length(spec2)) return(FALSE)
  return(all(sort(spec1) == sort(spec2)))
}

# Convert peptide (mass vector) to string representation
PeptideToString <- function(peptide) 
{
  return(paste(peptide, collapse = "-"))
}

# Main cyclopeptide sequencing algorithm
CyclopeptideSequencing <- function(spectrum) 
{
  candidate_peptides <- list(integer(0))  # Start with empty peptide
  final_peptides <- list()
  parent_mass <- ParentMass(spectrum)
  
  while (length(candidate_peptides) > 0) 
  {
    # Expand all candidate peptides
    candidate_peptides <- Expand(candidate_peptides)
    
    # Track peptides to remove
    to_remove <- c()
    
    for (i in seq_along(candidate_peptides)) 
    {
      peptide <- candidate_peptides[[i]]
      peptide_mass <- Mass(peptide)
      
      if (peptide_mass == parent_mass) 
      {
        # Check if cyclospectrum matches
        if (SpectraEqual(Cyclospectrum(peptide), spectrum)) 
        {
          peptide_str <- PeptideToString(peptide)
          
          # Check if not already in final peptides
          if (!any(sapply(final_peptides, function(p) p == peptide_str))) 
          {
            final_peptides[[length(final_peptides) + 1]] <- peptide_str
          }
        }
        to_remove <- c(to_remove, i)
      } else if (!IsConsistent(peptide, spectrum)) {
        to_remove <- c(to_remove, i)
      }
    }
    
    # Remove peptides that are complete or inconsistent
    if (length(to_remove) > 0) {
      candidate_peptides <- candidate_peptides[-to_remove]
    }
  }
  
  return(unlist(final_peptides))
}

# Example usage
# Spectrum from a known peptide (e.g., corresponding to "NQEL")
example_spectrum <- c(0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484)

result <- CyclopeptideSequencing(example_spectrum)
cat("Found peptides:\n")
print(result)

# Test with another example
cat("\n\nTesting with a simple spectrum:\n")
test_spectrum <- c(0, 71, 87, 158, 158, 229, 229, 316)
result2 <- CyclopeptideSequencing(test_spectrum)
print(result2)


# Speediness - BENCHMARKING:
###########################################################

# Benchmarking function for R
benchmark_algorithm <- function(spectrum, runs = 5) {
  times <- numeric(runs)
  
  for (i in 1:runs) {
    start_time <- Sys.time()
    result <- CyclopeptideSequencing(spectrum)
    end_time <- Sys.time()
    times[i] <- as.numeric(end_time - start_time, units = "secs")
  }
  
  list(
    result = result,
    avg_time = mean(times),
    min_time = min(times),
    max_time = max(times),
    times = times
  )
}

# Run benchmarks
cat("\n=== R Implementation Benchmarking ===\n\n")

# Example 1
cat("Example 1: NQEL spectrum\n")
example_spectrum <- c(0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484)
bench1 <- benchmark_algorithm(example_spectrum, runs = 5)
cat(sprintf("Average time: %.4f ms\n", bench1$avg_time * 1000))
cat(sprintf("Min time: %.4f ms\n", bench1$min_time * 1000))
cat(sprintf("Max time: %.4f ms\n\n", bench1$max_time * 1000))

# Example 2
cat("Example 2: Simple spectrum\n")
test_spectrum <- c(0, 71, 87, 158, 158, 229, 229, 316)
bench2 <- benchmark_algorithm(test_spectrum, runs = 5)
cat(sprintf("Average time: %.4f ms\n", bench2$avg_time * 1000))
cat(sprintf("Min time: %.4f ms\n", bench2$min_time * 1000))
cat(sprintf("Max time: %.4f ms\n\n", bench2$max_time * 1000))

# Example 3
cat("Example 3: Larger spectrum\n")
large_spectrum <- c(0, 71, 99, 101, 103, 128, 129, 170, 200, 202, 
                    228, 231, 257, 299, 303, 328, 330, 332, 333, 
                    403, 431, 434, 461, 534)
bench3 <- benchmark_algorithm(large_spectrum, runs = 3)
cat(sprintf("Average time: %.4f ms\n", bench3$avg_time * 1000))
cat(sprintf("Min time: %.4f ms\n", bench3$min_time * 1000))
cat(sprintf("Max time: %.4f ms\n", bench3$max_time * 1000))