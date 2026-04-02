#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 3 - Peptide Encoding problem - Highlight which DNA codons are encoding given peptide
# Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
# Output: All substrings of Text encoding Peptide (if any such substrings exist).
# NOTE: This is MUCH BETTER than Wk3_1. USE THIS for any future imports
#########################################################################################

# Constants
RNA_CODON_MAP = {
    "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

# Peptide encoding problem
def peptide_encoding(dna_string, peptide_string):
    """
    Takes a DNA string, uses both forward and reverse strands,
    transcribes to RNA  and finds codons that encode given peptide
    based upon the codon table (STOP is not translated)

    Parameters:
    dna_string -> Single stranded DNA string to test against
    peptide_string -> peptide string letters to find substring codons for

    Returns:
    List of substrings encoding peptide_string
    """
    answer = []
    peptide_len = len(peptide_string) * 3 # Because each peptide is 3 kmer long
    dna_len = len(dna_string)

    for i in range(dna_len - peptide_len + 1):
        substring = dna_string[i:i+peptide_len]

        # Forward strand
        rna_substring = dna_to_rna_transcription(substring)
        if rna_to_protein_translation(rna_substring) == peptide_string:
            answer.append(substring)

        # Reverse strand
        substring_rev = reverse_complement(substring)
        rna_substring_rev = dna_to_rna_transcription(substring_rev)
        if rna_to_protein_translation(rna_substring_rev) == peptide_string:
            answer.append(substring)

    return answer


# Helper functions
###########################
# 1: dna to rna map
def dna_to_rna_transcription(dna_string):
    """
    Takes a DNA string, and finds corresponding transcribed RNA

    Parameters:
    dna_string -> Single stranded DNA string to trnascribe

    Returns:
    RNA strand transcribed
    """
    #text = "AAAACCCGGT"
    rna_string = dna_string.replace("T", "U")
    return rna_string

# 2: reverse complement
def reverse_complement(dna_string):
    """
    Takes a DNA string, and finds reverse complement

    Parameters:
    dna_string -> Single stranded DNA string to reverse complement

    Returns:
    reverse complement
    """
    #text = "AAAACCCGGT"
    complements = {"A":"T", "T":"A", "C":"G", "G":"C"}
    complement = [complements[base] if base in complements else '' for base in dna_string]
    reverse_complement = ''.join(complement[::-1])
    return reverse_complement

# Protein translation
# String Reconstruction from patterns function
def rna_to_protein_translation(rna_string):
    """
    Takes an RNA string and tranlates it into a protein
    based upon the codon table (STOP is not translated)

    Parameters:
    rna_string -> RNA string to translate into protein letters

    Returns:
    String of protein letters
    """
    protein_string = ""

    # Loop through the RNA string in steps of 3 (codon length)
    for i in range(0, len(rna_string), 3):
        codon = rna_string[i:i+3]
        amino_acid = RNA_CODON_MAP.get(codon, '')

        # Ignore stop codon and codons less than 3 characters (prob not necessary to check)
        if amino_acid == 'STOP' or amino_acid == '' or len(codon) < 3:
            continue

        protein_string += amino_acid
    return protein_string


###########################################################################

if __name__ == "__main__":

# Sample test
    dna_string = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
    peptide_string = 'MA'
    # Expected answer = ATGGCC GGCCAT ATGGCC
    answer = peptide_encoding(dna_string, peptide_string)
    print('\n'.join(answer))

# # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
#     with open(file_path, 'r') as file:
#        dna_string = file.readline()
#        peptide_string = file.readline()

#     answer = peptide_encoding(dna_string, peptide_string)
#     # print('\n'.join(answer))

#     with open("Wk3_2_output.txt", "w") as output_file:
#         output_file.write('\n'.join(answer))


# # For Exercise Break

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: For Wk3_2 exercise break, we need to concatenate all the kmers on different lines
#     with open(file_path, 'r') as file:
#        dna_string = file.read().replace('\n', '')
#     # print(dna_string)

#     # Val-Lys-Leu-Phe-Pro-Trp-Phe-Asn-Gln-Tyr
#     peptide_string = 'VKLFPWFNQY'

#     answer = peptide_encoding(dna_string, peptide_string)
#     # print(len(answer))

#     with open("Wk3_2_exercisebreak_output.txt", "w") as output_file:
#         output_file.write('\n'.join(answer))
#         # Answer is 0!!