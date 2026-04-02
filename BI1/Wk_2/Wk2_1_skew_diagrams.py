##########################################################################################
# TASK:
# Create a skew diagram to identify ori. Count prior skew + 1 if G and prior skew -1 if C
# Fund the point (i.e., ori) where the value goes from decreasing to increasing OR
# Where the value is at a minimum
# Return the minimum value and the index
##########################################################################################


#CONSTANTS:
# TEXT = input("Please enter your text: ") #"CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
# file_name = input("Please enter file name: ")
# with open(file_name, 'r') as file:
#     GENOME = file.read()

# Original Genome
#GENOME = 'GAGCCACCGCGATA'
# Updated Genome
# GENOME = 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'

# REad file
file_name = input("Please enter file name: ")
with open(file_name, 'r') as file:
    GENOME = file.read()

# # Quiz:
# GENOME = 'GATACACTTCCCGAGTAGGTACTG' 

def min_skew(func):
    def wrapper(*args, **kwargs):
        func_answer = func(*args, **kwargs)
        #print(func_answer)
        # print('The skew is as follows ---->', ' '.join(map(str, func_answer)))
        min_skew_val = min(func_answer)
        min_skew_indices = [i for i, val in enumerate(func_answer) if val == min_skew_val]
        
        print("And the min skew is at nucleotide indices ----> ", ' '.join(map(str, min_skew_indices)))
        print("With the min value of ----> ", min_skew_val)
        # return min_skew_indices
    return wrapper

# BUILD SKEW
@min_skew
def skew(GENOME):
    skew_list = [0]
    current_skew = 0
    for i in GENOME:
        if i == 'G':
            current_skew += 1
        elif i == 'C':
            current_skew -= 1
        else:
            pass
        skew_list.append(current_skew)
    return skew_list

# PRETTY PRINT IT!

#print('The skew is as follows ---->', ' '.join(map(str, skew(GENOME))))


if __name__ == '__main__':
    skew(GENOME)
