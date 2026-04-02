# Step 1: Identify a kmer length pattern in a text string
# Step 2: Identify which of these pattern texts occurs the most frequently

# text = "TTGTTGTGTCAC"
# kmer_len = 2

text = input("Please enter your string: ")
kmer_len = int(input("Please enter the kmer length: "))

# Let'd develop a function to find the pattern of length k-mer and count the highest frequency pattern
def FrequencyTable(text, kmer_len):
    # Step 1: Identify a kmer length pattern in a text string
    freq_table = {}
    for i in range(len(text)-kmer_len+1):
        pattern = text[i:i+kmer_len]
        if pattern in freq_table:
            freq_table[pattern] += 1
        else:
            freq_table[pattern] = 1
    # Step 2: Identify which of these pattern texts occurs the most frequently
    max_count = max(freq_table.values())
    most_freq_patterns = [key for key, value in freq_table.items() if value == max_count]

    return freq_table, max_count, most_freq_patterns

results = FrequencyTable(text, kmer_len)
print(results)
print("Max frequency: ", results[1])
print("Patterns: ", results[2])