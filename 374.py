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
# def longestCommonSequence():  This function
def longestCommonSequence(sequence):
    pass


# Testing Area
# Main Function

