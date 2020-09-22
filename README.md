# Bioinformatics Homework 2 [CSCI 4800]

You are to implement a dynamic programming algorithm for pairwise sequence alignment: 
a program that does global alignment using a linear gap penalty function.

The program may need to operate on either DNA or protein sequences.

You can write your program in any language you like but the deliverable once compiled (or interpreted) should accept command line arguments. 
The Global alignment program should take the following as command-line arguments: 

    1 the name of a file containing the two sequences to be aligned, 
    2 the name of a file containing a substitution matrix, 
    3 an integer value for the gap penalty (g). 

We will ask that your program compute "highroad" alignments. 

For the global program, which uses a linear gap penalty function, you should determine this alignment as described in class. 

You should follow the convention that the rows of Dynamic Programming matrix correspond to the characters of the first sequence given, and the columns of the matrix correspond to characters of the second sequence given. 

## Input
You should run your global alignment program on the following pair of sequences

1. [human/mouse p53](https://drive.google.com/file/d/1nQ3MSouNVDC2aK5Xkk3vwJOoliEsl9Ru/view?usp=sharing)
   
Additionally, we will test your programs on several held-aside pairs of sequences.

For the global alignments that you turn in, you should set the gap penalty parameter g = -6.

For all of DNA sequence alignments, use [DNA substitution matrix](https://drive.google.com/file/d/1ytIw21nBuSiSqq5-gNJUgYEc0cdDHT_1/view?usp=sharing), and for protein sequence alignment runs, you should use the [Blosum-62](https://drive.google.com/file/d/1eghD2Q6Me8eDkfyCOKAENITer8aI-TtZ/view?usp=sharing) substitution matrix. 
You should assume that the sequence files given as input will be in the same format as those given above. That is, the first line of the file lists the first sequence, the second line is empty, the third line lists the second sequence, and subsequent lines can be ignored

## Output
Your programs should provide the output in the following format:
 
    First line: print the string “Global alignment”.

    Second line: Print “DNA Sequence” or “Protein Sequence”.

    Third line: Print the values of g

    Fourth, fifth and sixth lines: Print the pairwise sequence alignment result, where 4th and 6th lines denote the first and second aligned sequences, and the 5th line contain any matching symbol, like “|”. In the 4th and 6th lines, please use “-” character to denote gaps.

    Seventh line: Print the alignment score. 
    
## How to turn in Your HW
You should turn in the following items bundled in a ZIP archive:
1. Results
    * Human-mouse-alignment.txt
2. Source+Codes
    * *Your codes go here*
    * README.txt *[describing how to use your codes to get the above alignment results]*

