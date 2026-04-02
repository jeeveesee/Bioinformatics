# Get the reverse complement of a sample input - meeee did it!! 

# UNCOMMENT BELOW TO READ FROM FILE
#with open('Wk1_3_reversecomplements_dataset_exam2.txt') as reader:
#    text = reader.read()

##### print(text)

#text = "AAAACCCGGT"
complements = {"A":"T", "T":"A", "C":"G", "G":"C"}

# Read from input
text = input("Please enter your text: ")

def reverse_complement(text):
    answer = [complements[base] if base in complements else '' for base in text]
    return ''.join(answer[::-1])

result = reverse_complement(text)
print(result)
