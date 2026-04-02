#########################################################################################
# Genome Sequencing  - Bioinformatic III Course from Coursera
#
# Week 1 - Dynamic Programming Change creation
# Input: An integer money and an array Coins = (coin1, ..., coind)
# Output: The minimum number of coins with denominations Coins that changes money

"""
DPChange(money, Coins)
    MinNumCoins(0) ← 0
    for m ← 1 to money
        MinNumCoins(m) ← ∞
            for i ← 0 to |Coins| - 1
                if m ≥ coini
                    if MinNumCoins(m - coini) + 1 < MinNumCoins(m)
                        MinNumCoins(m) ← MinNumCoins(m - coini) + 1
    output MinNumCoins(money)
"""
# NOTE: YOU HAVE TO RUN THIS FROM THE BI3 FOLDER WITH THE COMMAND:
# python -m Wk1.Wk1_1_DPChange (NO .py)
#########################################################################################

def DPChange(money, Coins):
    """
    Uses dynamic programming to find the coins needed to make a given amount of money

    Parameters:
    money -> value of money that we need change for
    Coins -> list of coin denominations

    Returns:
    number of Coins needed to make change
    """
    global min_num_coins
    min_num_coins = [0]  # MinNumCoins(0) ← 0
    for m in range(1, money + 1):  # for m ← 1 to money
        min_num_coins.append(float('inf'))  # MinNumCoins(m) ← ∞
        # print(f"{min_num_coins=}")
        for coin in Coins:  # for i ← 0 to |Coins| - 1
            # print(f"{m=}")
            # print(f"{coin=}")
            if m >= coin:  # if m ≥ coini
                # print("m - coin =", m-coin)
                # print(f"min_num_coins[{m-coin}]=", min_num_coins[m-coin])
                # print(f"min_num_coins[{m-coin}] + 1 =", min_num_coins[m-coin] + 1)
                # print(f"min_num_coins[{m}] =", min_num_coins[m])
                if min_num_coins[m - coin] + 1 < min_num_coins[m]:  # if MinNumCoins(m - coini) + 1 < MinNumCoins(m)
                    min_num_coins[m] = min_num_coins[m - coin] + 1  # MinNumCoins(m) ← MinNumCoins(m - coini) + 1
                    # print(f"Yes, it fits, so min_num_coins[",{m},"] =", min_num_coins[m])
    return min_num_coins[money]  # output MinNumCoins(money)

###########################################################################

if __name__ == "__main__":

# # Sample test
#     money = 40
#     Coins_raw = '50 25 20 10 5 1'
#     Coins = list(map(int, Coins_raw.split()))
#     # Expected answer = 2
#     answer = DPChange(money, Coins)
#     # print('\n'.join(answer))
#     print(answer)

# For exam
    money = 25
    Coins_raw = '2 3'
    Coins = list(map(int, Coins_raw.split()))
    # Expected answer = 2
    answer = DPChange(money, Coins)
    # print('\n'.join(answer))
    print(answer)

# # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
#     with open(file_path, 'r') as file:
#         # Remove trailing newline/whitespace from the peptide string (was causing '\n' to be iterated)
#         money = file.readline().strip()
#         Coins_raw = file.readline()

#     Coins = list(map(int, Coins_raw.split()))

#     answer = DPChange(int(money), Coins)
#     # print('\n'.join(answer))
#     print(answer)

#     with open("Wk1_1_output.txt", "w") as output_file:
#         output_file.write(str(answer))


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

# # EXAM
#     peptide_string = 'PEEP'
#     experimental_spectrum = [0, 97, 129, 129, 129, 194, 226, 323, 323, 355, 452]
#     answer = linear_scoring(peptide_string, experimental_spectrum)
#     # print('\n'.join(answer))
#     print(answer)