# Identify a pattern in a text string and then
# identify which text occurs the most frequently

text = "TTGTGAACCTTGTCCGTTGTGAACCTTGTCCGGCTGTACTTTGTGAACTCCGGTACTTGTCCGGTTAAAATCCGGTATTGTGAACGTTAAAAGCTGTACTGTTAAAAGTTAAAAGTTAAAAGCTGTACTCTTGTCCGTTGTGAACTCCGGTATCCGGTATTGTGAACCTTGTCCGGCTGTACTGTTAAAACTTGTCCGTTGTGAACTTGTGAACTCCGGTAGTTAAAACTTGTCCGTCCGGTAGCTGTACTCTTGTCCGTCCGGTAGCTGTACTGCTGTACTGTTAAAAGCTGTACTGCTGTACTCTTGTCCGGTTAAAAGTTAAAAGTTAAAAGTTAAAATCCGGTAGCTGTACTGCTGTACTTTGTGAACGTTAAAATCCGGTAGTTAAAATTGTGAACTTGTGAACTTGTGAACGTTAAAAGCTGTACTGTTAAAATTGTGAACTTGTGAACTTGTGAACGTTAAAATCCGGTAGTTAAAACTTGTCCGGCTGTACTGCTGTACTTTGTGAACGCTGTACTTTGTGAACTTGTGAACGCTGTACTTCCGGTAGTTAAAATCCGGTAGCTGTACTTTGTGAACCTTGTCCGGTTAAAATCCGGTATCCGGTACTTGTCCGGCTGTACTGCTGTACTGTTAAAACTTGTCCGCTTGTCCGTCCGGTATCCGGTACTTGTCCGTTGTGAACCTTGTCCGGCTGTACTGTTAAAAGTTAAAACTTGTCCGGCTGTACTCTTGTCCGCTTGTCCGTCCGGTATTGTGAACGCTGTACTGTTAAAATTGTGAACTCCGGTACTTGTCCG"
kmer_len = 14

# Count the number of a single k-mer patterns
def PatternCount(text, pattern):
    count = 0
    for i in range(len(text) - len(pattern)+1):
        if text[i:i+len(pattern)] == pattern:
            count += 1
    return count

# Count all k-mers
pattern_list = []
count_list = []
def FrequentWords(text, kmer_len):
    for i in range(len(text)-kmer_len):
        pattern = text[i:i+kmer_len]
        pattern_list.append(pattern)
        count_list.append(PatternCount(text, pattern))
    max_count = max(count_list)
    max_index = [index for index, value in enumerate(count_list) if value == max_count]
    return max_count, max_index

results = FrequentWords(text, kmer_len)
# print(results)
# print(pattern_list)
print("Max frequency: ", results[0])
print("Max value indices: ", results[1])
kmers = set(pattern_list[i] for i in results[1])
print("Patterns:", kmers)
print("Patterns broken by space :-) :", *kmers)

print(results)
# print("max count is:", results.max_count)
# print("index for max counts are:", results.max_index)
