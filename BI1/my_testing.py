


alpha = []
x = {'a': 2, 'b': 4, 'c':  5, 'd': 1}
y = {'e': 2, 'f': 4, 'g':  5, 'h': 1}
z = {'i': 2, 'j': 4, 'k':  5, 'l': 1}
t = 3

listula = [x, y, z]

for tesla in listula:
    [alpha.append(kmer) for kmer in tesla.keys() if tesla[kmer] >= t]

print(alpha)




# for i in range(10):
#     print(i)
#     z = 
#     [x.append(y) for y in ]



# def find_clumps(text, kmer_len, window_len, kmer_freq):
#     kmer_clumps = []
#     text_len = len(text)

#     for i in range(text_len - window_len):
#         window = text[i:window_len]
#         freqmap = frequency_table(window, kmer_len)[0] # Just return the freq_table
#         # TODO: In the below, we should try to capture the length as well
        
        
#         [kmer_clumps.append(kmer) for kmer in freqmap.keys() if freqmap[kmer] >= kmer_freq]
    
#     return set(kmer_clumps)
