#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 3 - RNA to Protein Translation
# Input: An RNA string Pattern and the array GeneticCode.
# Output: The translation of Pattern into an amino acid string Peptide.
#########################################################################################

# String Reconstruction from patterns function
def rna_to_protein_translation(rna_codon_map, rna_string):
    """
    Takes an RNA string and tranlates it into a protein
    based upon the codon table (STOP is not translated)

    Parameters:
    rna_codon_map -> RNA Codon (in Python dictionary format)
    rna_string -> RNA string to translate into protein letters

    Returns:
    String of protein letters
    """
    protein_string = ""

    # Loop through the RNA string in steps of 3 (codon length)
    for i in range(0, len(rna_string), 3):
        codon = rna_string[i:i+3]
        amino_acid = rna_codon_map.get(codon, '')

        # Ignore stop codon and codons less than 3 characters (prob not necessary to check)
        if amino_acid == 'STOP' or amino_acid == '' or len(codon) < 3:
            continue

        protein_string += amino_acid
    return protein_string


############################################################################

if __name__ == "__main__":

    rna_codon_map = {
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

# Sample test
    rna_string = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
    # Expected answer = MAMAPRTEINSTRING
    # rna_string = "CCAAGUACAGAGAUUAAC"

    answer = rna_to_protein_translation(rna_codon_map, rna_string)
    print(answer)

# # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     with open(file_path, 'r') as file:
#        rna_string = file.read()

#     answer = rna_to_protein_translation(rna_codon_map, rna_string)
#     # print(answer)

#     with open("Wk3_1_output.txt", "w") as output_file:
#         output_file.write(answer)