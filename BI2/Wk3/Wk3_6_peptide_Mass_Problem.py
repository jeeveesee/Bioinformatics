#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 3 - peptides with given mass problem
# Input: An integer m.
# Output: The number of linear peptides having integer mass m.

# This requires dynamic programming
# dp[0] = 1 means that there is 1 way to make a mass of 0 i.e., using no AA
# dp[mass] = dp[mass] + dp[mass-mass of amino acid]
# Essentially, to get to mass, you need to see how many sequences help you
# get to mass -a... e.g., if the last aa has a mass of 5, you see how many
# sequences can get you to mass-5, then add how many sequences get you to mass
# then add them up (for example, to get to mass, you might just be able to 
# get it with one aa, which is the same as the total mass!)

#########################################################################################

AMINO_ACID_MASSES = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

def count_peptides_with_mass(target_mass):
	dp = [0] * (target_mass + 1)
	dp[0] = 1  # base case: one way to make mass 0 (empty peptide)

	for mass in range(1, target_mass + 1):
		#print(f"\n{mass=}")
		for a in AMINO_ACID_MASSES:
			#print(f"{a=}")
			#print(f"{mass-a=}")
			if mass - a >= 0:
				#print(f"{dp[mass]=}")
				#print(f"{dp[mass-a]=}")
				dp[mass] += dp[mass - a]
				#print(f"New {dp[mass]=}")
	return dp[target_mass]

###########################################################################

if __name__ == "__main__":

# # Sample test
#     target_mass = 1024
#     # Expected answer = 14712706211
#     answer = count_peptides_with_mass(target_mass)
#     ## print(''.join(answer))
#     print(answer)


# # From file

    # Get dataset
    from pathlib import Path as partho
    current_dir = partho(__file__).parent
    filename = input("Please enter the filename: ")
    file_path = current_dir / filename

    #NOTE: Make sure to remove the extra line at the end of the file
    with open(file_path, 'r') as file:
       target_mass = file.readline()

    answer = count_peptides_with_mass(int(target_mass))
    # print('\n'.join(answer))

    with open("Wk3_6_output.txt", "w") as output_file:
    #   output_file.write('\n'.join(answer))
    #   output_file.write(" ".join(map(str, answer)))
        output_file.write(str(answer))



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

#     answer = peptide_encoding(RNA_CODON_MAP, dna_string, peptide_string)
#     # print(len(answer))

#     with open("Wk3_2_exercisebreak_output.txt", "w") as output_file:
#         output_file.write('\n'.join(answer))
#         # Answer is 0!!