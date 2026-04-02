#########################################################################################
# Genome Sequencing  - Bioinformatic II Course from Coursera
#
# Week 4 - Convolution Cyclopeptide Sequencing with real Tyrocidine B1 spectrum

"""
Try to sequence the tyrocidine corresponding to the real experimental spectrum below.
Since the fragmentation technology used for generating the spectrum tends to produce ions with charge +1,
you can safely assume that all charges are +1. Return the peptide as a collection of space-separated integer masses.
"""

# NOTE: YOU HAVE TO RUN THIS FROM THE BI2 FOLDER WITH THE COMMAND:
# python -m Wk4.Wk4_9_Real_TyrocidineB1_Sequencing (NO .py)
#########################################################################################

# Imports:
from Wk4.Wk4_8_Convolution_Cyclopeptide_Sequencing import convolution_cyclopeptide_sequencing

# def tyrocidineB1_sequencing(experimental_spectrum, M, N):
#     leaders, best_score = convolution_cyclopeptide_sequencing(experimental_spectrum, M, N)
#     return leaders, best_score




###########################################################################

if __name__ == "__main__":

# # Sample test
#     experimental_spectrum_raw = "57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493"
#     experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))
#     M = 20
#     N = 60
#     # Expected answer = 99-71-137-57-72-57
#     leaders, best_score = convolution_cyclopeptide_sequencing(experimental_spectrum, M, N)
#     formatted = ["-".join(map(str, p)) for p in leaders]
#     print("Leader peptides are:\n")
#     print("\n".join(formatted))
#     print("\nBest score is:", best_score)
#     # print(*answer)

# # From file

#     # Get dataset
#     from pathlib import Path as partho
#     current_dir = partho(__file__).parent
#     filename = input("Please enter the filename: ")
#     file_path = current_dir / filename

#     #NOTE: For Wk3_2, there is a newline character after the peptide string, so remove that
#     with open(file_path, 'r') as file:
#         # Read peptide string and remove ALL whitespace (spaces/newlines) — dataset may contain spaces
#         # leaderboard = file.readline().split()
#         M = int(file.readline())
#         N = int(file.readline())
#         experimental_spectrum_raw = file.readline()
#     experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))

#     leaders, best_score = convolution_cyclopeptide_sequencing(experimental_spectrum, M, N)
#     formatted = ["-".join(map(str, p)) for p in leaders]
#     answer = "\n".join(formatted)

#     with open("Wk4/Wk4_8_output.txt", "w") as output_file:
#         # output_file.write(str(answer))
#         output_file.write(answer)


# # For Exercise Break
    # For Wk4_9  Tyrocidine B1 Spectrum
    M = 20
    N = 1000
    experimental_spectrum_raw = "371.5 375.4 390.4 392.2 409.0 420.2 427.2 443.3 446.4 461.3 471.4 477.4 491.3 505.3 506.4 519.2 536.1 546.5 553.3 562.3 588.2 600.3 616.2 617.4 618.3 633.4 634.4 636.2 651.5 652.4 702.5 703.4 712.5 718.3 721.0 730.3 749.4 762.6 763.4 764.4 779.6 780.4 781.4 782.4 797.3 862.4 876.4 877.4 878.6 879.4 893.4 894.4 895.4 896.5 927.4 944.4 975.5 976.5 977.4 979.4 1005.5 1007.5 1022.5 1023.7 1024.5 1039.5 1040.3 1042.5 1043.4 1057.5 1119.6 1120.6 1137.6 1138.6 1139.5 1156.5 1157.6 1168.6 1171.6 1185.4 1220.6 1222.5 1223.6 1239.6 1240.6 1250.5 1256.5 1266.5 1267.5 1268.6"
    # experimental_spectrum = list(map(int, experimental_spectrum_raw.split()))
    # Given the +1 charge, the mass/charge ratio i.e., m/z = (mass of fragment + mass of proton)/charge
    # As such we need to reduce the mass by 1 before rounding
    experimental_spectrum_raw_floats = [float(x) for x in experimental_spectrum_raw.split()]
    experimental_spectrum = [int(round(f - 1.0)) for f in experimental_spectrum_raw_floats]
    print(experimental_spectrum)
    leaders, best_score = convolution_cyclopeptide_sequencing(experimental_spectrum, M, N)
    formatted = [" ".join(map(str, p)) for p in leaders]
    answer = "\n".join(formatted)
    print(answer)