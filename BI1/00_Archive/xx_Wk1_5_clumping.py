# Identify a pattern in a text string and then identify the indices where it matches

text = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
k = 5
l = 50
t = 4

# Count the number of single k-mer patterns in a given text
def PatternCount(text_to_check, pattern, t):
    count = 0
    for i in range(len(text_to_check) - len(pattern)+1):
        if text_to_check[i:i+len(pattern)] == pattern:
            count += 1
    return count if count >= t else None

# Create a window, slide window across the text to count occurences
# Then move the window by one to repeat above again with same pattern
# Then restart with the window moved by one and pick up new pattern and repeat
text_len = len(text)


def clumping(text, k, l, t):
    for i in range(len(text) - k):
        pattern = text[i:i+k]
        text_len = len(text)
        for j in range(min(text_len, l)):
            text_to_check = text[j:j+min(text_len, l)+1]
            text_len -= 1
            print('Pattern', pattern, ' Index', i+j, ' Count', PatternCount(text_to_check, pattern, t))


clumping(text, k, l, t)
#
#
#
#
#
# # Count all k-mers
# pattern_list = []
# count_list = []
# def FrequentWords(text, k, l, t):
#     for i in range(len(text)-min(l, len(text))):
#         pattern = text[i:i+k]
#         pattern_list.append(pattern)
#         for j in range(min(l, len(text))):
#             count_list.append(PatternCount(text, l, pattern))
#     max_index = [index for index, value in enumerate(count_list) if value == t]
#     return t, max_index
#
# results = FrequentWords(text, k, l, t)
# print(results)
# print(pattern_list)
# print("Max frequency: ", results[0])
# print("Max value indices: ", results[1])
# kmers = set(pattern_list[i] for i in results[1])
# print("Patterns:", kmers)
# print("Patterns broken by space :-) :", *kmers)
