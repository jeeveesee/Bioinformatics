##########################################################################################
# TASK:
# Identify a pattern in a text string and then identify how many times that shows up in a given window
##########################################################################################


#CONSTANTS:
# TEXT = input("Please enter your text: ") #"CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
file_name = input("Please enter file name: ")
with open(file_name, 'r') as file:
    TEXT = file.read()
KMER_LEN = 9
WINDOW_LEN =500
KMER_FREQ = 3

##########################################################################################
# STEP 1:
# Count the number of single k-mer patterns in a given text
# This is our frequencytable function from Wk1_2
##########################################################################################
def frequency_table(text, kmer_len):
    # Step 1: Identify a kmer length pattern in a text string
    freq_table = {}
    for i in range(len(text)-kmer_len+1):
        pattern = text[i:i+kmer_len]
        if pattern in freq_table:
            freq_table[pattern] += 1
        else:
            freq_table[pattern] = 1
    # Step 2: Identify which of these pattern texts occurs the most frequently
    # DON'T REALLY NEED THESE FOR THIS CLUMPING PROBLEM. COMMENTING FOR EFFICIENCY
    # max_count = max(freq_table.values())
    # most_freq_patterns = [key for key, value in freq_table.items() if value == max_count]

    return freq_table#, max_count, most_freq_patterns



##########################################################################################
# STEP 2:
# Create a window, slide window across a specific length of the big text (e.g., 20 BPs window 
# over 500 big text length). # to build frequency table map. Then move the window by one to 
# repeat above again with same pattern. during each run, count if a particular kmer was found 
# more than t times in that length"""
##########################################################################################

def find_clumps(text, kmer_len, window_len, kmer_freq):
    kmer_clumps = []
    text_len = len(text)

    for i in range(text_len - window_len):
        window = text[i:i+window_len]
        freqmap = frequency_table(window, kmer_len)#[0] # Just return the freq_table
        # TODO: In the below, we should try to capture the length as well
        [kmer_clumps.append(kmer) for kmer in freqmap.keys() if freqmap[kmer] >= kmer_freq]
    
    return set(kmer_clumps)



##########################################################################################
# TEST IT:
##########################################################################################

# kmer_clumps = []
# pre_answer = frequency_table(TEXT, KMER_LEN)[0]
# print(pre_answer)
# answer = [kmer for kmer in pre_answer.keys() if pre_answer[kmer] >= KMER_FREQ]
# print(answer)

# print(find_clumps(TEXT, KMER_LEN, WINDOW_LEN, KMER_FREQ))
answer = find_clumps(TEXT, KMER_LEN, WINDOW_LEN, KMER_FREQ)
print(len(answer))
print(answer)
# print(frequency_table(TEXT, KMER_LEN)[2])
