'''
CSE 374 Group Assignment
Names : Alex Habegger,

PLEASE INPUT YOUR NAME SO THE GROUP KNOWS THAT YOU VIEWED THE CODE
PLEASE INPUT YOUR NAME SO THE GROUP KNOWS THAT YOU VIEWED THE CODE
PLEASE INPUT YOUR NAME SO THE GROUP KNOWS THAT YOU VIEWED THE CODE

Rough Code With Documentation
'''

# ALGORITHM ONE
# def translation():  This function takes DNA sequence as the only
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

   F1 = translation(sequence)
   F2 = translation(sequence[1:])
   F3 = translation(sequence[2:])
   F4 = translation(rna)
   F5 = translation(rna[1:])
   F6 = translation(rna[2:])
   return (F1, F2, F3, F4, F5, F6)



# ALGORITHM THREE
# ALGORITHM THREE
# def longestCommonSequence(): This function acts as the
# main calling fucntion
def longestCommonSequenceMain(seq1, seq2):
    return LCSubStr(seq1, seq2, len(seq1), len(seq2))


# Python3 implementation of Finding
# Length of Longest Common Substring

# Returns length of longest common
# substring of X[0..m-1] and Y[0..n-1]


def LCSubStr(X, Y, m, n):
    # Create a table to store lengths of
    # longest common suffixes of substrings.
    # Note that LCSuff[i][j] contains the
    # length of longest common suffix of
    # X[0...i-1] and Y[0...j-1]. The first
    # row and first column entries have no
    # logical meaning, they are used only
    # for simplicity of the program.

    # LCSuff is the table with zero
    # value initially in each cell
    LCSuff = [[0 for k in range(n + 1)] for l in range(m + 1)]

    # To store the length of
    # longest common substring
    result = 0

    # Following steps to build
    # LCSuff[m+1][n+1] in bottom up fashion
    for i in range(m + 1):
        for j in range(n + 1):
            if (i == 0 or j == 0):
                LCSuff[i][j] = 0
            elif (X[i - 1] == Y[j - 1]):
                LCSuff[i][j] = LCSuff[i - 1][j - 1] + 1
                result = max(result, LCSuff[i][j])
            else:
                LCSuff[i][j] = 0
    return result



# ALGORITHM FOUR
def LongestGeneExpression(frames1, frames2):
    print("Longest Gene Expression")
    framePair = ["",""]
    longestLength = 0
    for frame1 in frames1:
        for frame2 in frames2:
            currentLength = longestCommonSequenceMain(frame1, frame2)
            if (currentLength) > (longestLength):
                longestLength = currentLength
                framePair[0] = frame1
                framePair[1] = frame2

    print(longestLength)
    print(framePair[0])
    print(framePair[1])


# Testing Area
# Main Function
if __name__ == "__main__":
    # Example
    dna_seq1 = "AATTGGGGATCGATCGCATCAGCTAGCATCGACTAGCTAGC"
    dna_seq2 = "AGTTGGGGCCCGTTGTGCAAAGTTCGCTAATCGACTAGCTAGC"


    frames1 = sixFrame(dna_seq1)
    frames2 = sixFrame(dna_seq2)

    print("DNA ONE ", dna_seq1)
    print("Six Frames")
    for frame in frames1:
        print(frame)

    print("\n\nDNA TWO ", dna_seq2)
    print("Six Frames")
    for frame in frames2:
        print(frame)

    print()
    LongestGeneExpression(frames1, frames2)
