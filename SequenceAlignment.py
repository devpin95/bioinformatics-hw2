#!/usr/bin/python

import argparse

SCORE = 'score'
DIR = 'dir'
DIAG = 'diagonal'
LEFT = 'left'
UP = 'up'

s1 = ''
s2 = ''
nmmatrix = []
submatrix = {}
gap = None


def vprint(s, v):
    if v:
        print(s)


def match_chars(c1, c2):
    value = submatrix[c1][c2]
    match = False
    if c1 == c2:
        match = True
    return value, match


def cost(i, j):
    print(i-1, j-1)
    s1_char = s1[i-1]
    s2_char = s2[j-1]

    # check if diag has been found
    # print(nmmatrix[i-1][j-1], nmmatrix[i-1][j], nmmatrix[i][j-1])

    if nmmatrix[i-1][j-1][SCORE] is None:
        cost(i-1, j-1)

    # check if up has been found
    if nmmatrix[i-1][j][SCORE] is None:
        cost(i - 1, j)

    # check if up has been found
    if nmmatrix[i][j-1][SCORE] is None:
        cost(i, j-1)

    diag = nmmatrix[i - 1][j - 1] + match_chars(s1_char, s2_char)
    smax = {'value': diag, 'dir': DIAG}

    leftv = nmmatrix[i][j - 1] + gap
    if leftv > smax['value']:
        smax['value'] = leftv
        smax['dir'] = LEFT

    upv = nmmatrix[i-1][j] + gap
    if upv > smax['value']:
        smax['value'] = upv
        smax['dir'] = UP

    nmmatrix[i][j][SCORE] = smax['value']
    nmmatrix[i][j][DIR] = smax['dir']


def do_align(s, m, g, v):
    print('Global alignment')  # required print
    # vprint("Sequence File: " + s, v)
    # vprint("Matrix File: " + m, v)
    # vprint("Gap Penalty: " + str(g), v)
    # vprint("Verbose: " + str(v), v)

    gap = g

    vprint("|Preparing to do alignment...", v)

    sfile = open(s, 'r')
    is_protein_seq = False

    if sfile:
        vprint('|\t' + s + " opened", v)

        s1 = sfile.readline().rstrip()  # read in the first sequence on line 1
        s1_n = len(s1)
        sfile.readline()  # read in the blank line on line 2
        s2 = sfile.readline().rstrip()  # read in the second sequence on line 3
        s2_n = len(s2)

        # Print the first 50 chars of both sequences
        vprint("|\t\tSequence 1 (" + str(s1_n) + "): " + s1[: 50] + "...", v)
        vprint("|\t\tSequence 2 (" + str(s2_n) + "): " + s2[: 50] + "...", v)

    else:
        print("Failed to open " + s, v)
        return

    vprint("|\tDetermining sequence type (DNA/Protein)...", v)

    if 'm' or 'M' or 'r' or 'R' or 'w' or 'W' in s1:
        is_protein_seq = True

    vprint("|\t\t-> Protein Sequence" if is_protein_seq else "|\t\t-> DNA Sequence", v)

    print("Protein Sequence" if is_protein_seq else "DNA Sequence")  # required print
    print(g)  # required print

    vprint("|\tReading matrix...", v)

    mfile = open(m, 'r')

    if mfile:
        vprint('|\t\tMatrix file ' + m + " opened", v)
        vprint('|\t\t\tConstructing substitution matix...', v)

        alphabet = mfile.readline().rstrip().split()  # read the top line of characters in the matrix
        alphabet = [char.lower() for char in alphabet]  # convert chars to lowercase
        alphabet_size = len(alphabet)

        vprint("|\t\t\tAlphabet (" + str(alphabet_size) + ") " + str(alphabet), v)

        for char in alphabet:
            submatrix[char.lower()] = {}

        for i in range(0, alphabet_size):
            isubvals = mfile.readline().rstrip().split()
            columnchar = isubvals[0].lower()  # remove the char from the front of the line, use it as first entry
            isubvals.pop(0)

            for j in range(0, len(isubvals)):
                rowchar = alphabet[j].lower()
                submatrix[columnchar][rowchar] = int(isubvals[j])

        vprint('|\t\tSubstitution matrix created', v)
        vprint('|\t\t\t-> ' + str(submatrix), v)

    else:
        print("Failed to open " + m)
        return

    vprint('|\n|\tCreating Needleman-Wuncsh matrix...', v)
    # vprint('|\t\t' + str(len(s1) + 1) + "x" + str(len(s2) + 1), v)

    for i in range(0, len(s1)+1):
        nmmatrix.append([])

    for li in nmmatrix:
        for i in range(0, len(s2)+1):
            tup = {'score': None, 'dir': None}
            li.append(tup)

    nmmatrix[0][0][SCORE] = 0  # Initialize the first number

    # print(nmmatrix)
    vprint('|\t\t' + str(len(nmmatrix)) + "x" + str(len(nmmatrix[0])), v)

    vprint('|\t\tInitializing row one with gap penalty...', v)
    for i in range(1, s2_n+1):
        nmmatrix[0][i][SCORE] = nmmatrix[0][i-1][SCORE] + g
        nmmatrix[0][i][DIR] = UP
    vprint('|\t\t\t-> ' + str(nmmatrix[0]) + '...', v)

    vprint('|\t\tInitializing column one with gap penalty...', v)
    for i in range(1, s1_n + 1):
        nmmatrix[i][0][SCORE] = nmmatrix[i-1][0][SCORE] + g
        nmmatrix[i][0][DIR] = LEFT

    vprint('|\n|\tCalculating costs...', v)

    cost(s1_n, s2_n)
    print(nmmatrix[s1_n+1][s2_n+1])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('sequences', type=str,
                        help='a file containing two sequences')
    parser.add_argument('matrix', type=str,
                        help='a file containing DNA substitution matrix or Blosum-62 substitution matrix')
    parser.add_argument('gpenalty', type=int,
                        help='gap penalty')
    parser.add_argument('--verbose', action="store_true",
                        help='Flag to print intermediate steps of algorithm')

    args = parser.parse_args()

    sequences = args.sequences
    matrix = args.matrix
    gap = args.gpenalty
    verbose = args.verbose
    do_align(sequences, matrix, gap, verbose)
