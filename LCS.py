'''
CSE 374 Group Assignment
Names : Alex Habegger, Josh Lawson, Michael Glum
Names : Dillon Watkins, Max Zaremba
'''

# Importing time in order to do time testing
import time

# Importing random in order to create random DNA sequences
import random


# ALGORITHM ONE
# def translation(sequence): This function takes DNA sequence as the only
# string argument. It will return a protein sequence as the returned
# string. The Dictionary is hard coded DNA sequences to Protein.
def translation(sequence):
    # Hard Coded Dictionary of Codons to Protiens
    # Reference : https://www.chemguide.co.uk/organicprops/aminoacids/dna5.html
    trans_dic = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

    rna_seq = sequence.replace('T', 'U')
    protiens = []

    # Counts by 3 to grab all codons in sequence
    for x in range(0, len(rna_seq) - 2, 3):
        # Some Sequence have N standing for unknown nucleotide
        # We will skip over these codons
        if 'N' not in rna_seq[x:x+3]:
            protiens.append(trans_dic[rna_seq[x:x+3]])

    return "".join(protiens)


# ALGORITHM TWO
# def sixFrame():  This function takes DNA sequence as the only string
# argument. It will return 6 translation frames. It will invoke
# translation() 6 times to conduct the 6-frame translation.
def sixFrame(sequence):
   rna = sequence.replace("T","a").replace("A","t").replace("G","c").replace("C","g")
   rna = rna[::-1].upper()

   # Creates the translations of possible binding sites
   F1 = translation(sequence)
   F2 = translation(sequence[1:])
   F3 = translation(sequence[2:])
   F4 = translation(rna)
   F5 = translation(rna[1:])
   F6 = translation(rna[2:])

   # Returns the six frames that will be used for analysis
   return (F1, F2, F3, F4, F5, F6)


# ALGORITHM THREE
# def longestCommonSequence():
# Modified version of Longest Common Substring but
# will not return a LCS with end codon '*'
def longestCommonSubstring(seq1, seq2):
    # Creates a table in order to be used during
    # Dynamic programming using a Tabulation method.
    table = [[0 for a in range(len(seq2) + 1)] for b in range(len(seq1) + 1)]

    # Creates result to store the length
    # of the current Longest Common Substring and
    # creates variables to store index of cell in
    # table.
    result = 0
    row = 0
    col = 0

    # Following steps to build
    for i in range(len(seq1) + 1):
        for j in range(len(seq2) + 1):
            if (i == 0 or j == 0):
                table[i][j] = 0
            elif (seq1[i - 1] == seq2[j - 1]):
                # Will not return a LCS with end codon '*'
                if (seq1[i - 1] != '*'):
                    table[i][j] = table[i - 1][j - 1] + 1
                    if result < table[i][j]:
                        result = table[i][j]
                        row = i
                        col = j
            else:
                table[i][j] = 0


    # allocate space for the LCS
    LCS = ['0'] * result

    while table[row][col] != 0:
        result -= 1
        LCS[result] = seq1[row - 1]

        row -= 1
        col -= 1

    # required longest common substring
    return ''.join(LCS)


# ALGORITHM FOUR
# Compares every frame in frame one to every frame in frame two
# to find the longest gene expression of all frame combinations
def LongestGeneExpression(frames1, frames2):
    framePair = ["",""]
    longestLength = 0
    longestGeneExpression = ''

    for frame1 in frames1:
        for frame2 in frames2:
            lcs = longestCommonSubstring(frame1, frame2)
            currentLength = len(lcs)
            if (currentLength) > (longestLength):
                longestLength = currentLength
                framePair[0] = frame1
                framePair[1] = frame2
                longestGeneExpression = lcs

    print("Longest Gene Expression Results")
    print("LGE ({1}) : {0}".format(longestGeneExpression, longestLength))
    print("Frame 1.{2} ({1}) : {0}".format(framePair[0], len(framePair[0]), frames1.index(framePair[0]) + 1))
    print("Frame 2.{2} ({1}) : {0}".format(framePair[1], len(framePair[1]), frames2.index(framePair[1]) + 1))


# Purely for Testing Helper
# Generates a random nucleotide sequence
def generateRandom(length):
    sequence = ''
    for x in range(length):
        nuc = ''
        # Chooses a random numbeer between 1 and 4
        num = random.randint(1, 4)
        if num == 1:
            nuc = 'A'
        if num == 2:
            nuc = 'T'
        if num == 3:
            nuc = 'C'
        if num == 4:
            nuc = 'G'
        sequence += nuc

    # Return the random nucleotide sequence
    return sequence


# Main Function
# Prints details of results of other Algos
def main(dna_seq1, dna_seq2):
    # Records start time
    time1 = time.time()

    # Creates six frames for two DNA sequences
    frames1 = sixFrame(dna_seq1)
    frames2 = sixFrame(dna_seq2)

    # Calls Algorithm Four with the two sets of frames as parameters
    LongestGeneExpression(frames1, frames2)

    # Records end time and prints elapsed time
    time2 = time.time()
    print("Time Elapsed : {:.5f} seconds\n".format(time2 - time1))


# Testing Area
if __name__ == "__main__":
    random.seed(7124)

    # FULL TESTING SET
    print(100)
    dna_seq1 = generateRandom(100)
    dna_seq2 = generateRandom(100)
    main(dna_seq1, dna_seq2)

    print(250)
    dna_seq1 = generateRandom(250)
    dna_seq2 = generateRandom(250)
    main(dna_seq1, dna_seq2)

    print(500)
    dna_seq1 = generateRandom(500)
    dna_seq2 = generateRandom(500)
    main(dna_seq1, dna_seq2)

    print(750)
    dna_seq1 = generateRandom(750)
    dna_seq2 = generateRandom(750)
    main(dna_seq1, dna_seq2)

    print(1000)
    dna_seq1 = generateRandom(1000)
    dna_seq2 = generateRandom(1000)
    main(dna_seq1, dna_seq2)

    print(2500)
    dna_seq1 = generateRandom(2500)
    dna_seq2 = generateRandom(2500)
    main(dna_seq1, dna_seq2)

    print(5000)
    dna_seq1 = generateRandom(5000)
    dna_seq2 = generateRandom(5000)
    main(dna_seq1, dna_seq2)

    print(7500)
    dna_seq1 = generateRandom(7500)
    dna_seq2 = generateRandom(7500)
    main(dna_seq1, dna_seq2)

    print(10000)
    dna_seq1 = generateRandom(10000)
    dna_seq2 = generateRandom(10000)
    main(dna_seq1, dna_seq2)
