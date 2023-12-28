# CSE 374 Group Assignment
### Alexander J. Habegger

This repository contains a Python implementation of algorithms to find the Longest Gene Expression (LGE) between two DNA sequences. The group members are Alex Habegger, Josh Lawson, Michael Glum, Dillon Watkins, and Max Zaremba.

## Features
* Random DNA sequence generator
* Translation of DNA sequences to protein sequences
* 6-frame translation of DNA sequences
* Longest Common Substring with modified end codon
* Longest Gene Expression (LGE) comparison between two sets of translated frames

## Usage
1. Clone the repository
2. Navigate to the directory containing the script
3. Run the script using python3 <script_name.py>
4. 
## Algorithms
1. translation(sequence) - Translates a given DNA sequence into a protein sequence using a hard-coded dictionary of codons to proteins.
2. sixFrame(sequence) - Takes a DNA sequence as input and returns six translation frames.
3. longestCommonSubstring(seq1, seq2) - Finds the Longest Common Substring (LCS) between two protein sequences, with a modification to not return an LCS with an end codon '*'.
4. LongestGeneExpression(frames1, frames2) - Compares every frame in the first set of frames to every frame in the second set of frames to find the Longest Gene Expression (LGE) of all frame combinations.

## Testing
The script contains a main function that generates random DNA sequences of varying lengths and then applies the implemented algorithms to find the Longest Gene Expression (LGE) between the generated sequences. The elapsed time for each test is also printed.

## Disclaimer
This project is part of an academic exercise and is not intended for commercial use.
