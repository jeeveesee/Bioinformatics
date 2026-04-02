# Gte the reverse complement of a sample input

with open('Wk1_3_reversecomplements_dataset.txt') as reader:
    text = reader.read()
# print(text)

#text = "AAAACCCGGT"
complements = {"A":"T", "T":"A", "C":"G", "G":"C"}

def reverse_complement(text):
    answer = ''
    for i in text:
        answer += complements[i]
    return answer[::-1]

result = reverse_complement(text)
print(result)
