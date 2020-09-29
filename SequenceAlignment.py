#!/usr/bin/python

import argparse

SCORE = 'score'
DIR = 'dir'
MATCH = 'match'
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


def traceback(i, j):
    ts = ''  # top sequence
    bs = ''  # bottom sequence
    ms = ''  # match sequence (bars to show a matching value)

    while j > 0 or i > 0:
        s1_char = '-'
        s2_char = '-'
        if j > 0:
            s2_char = s2[j - 1]
        if i > 0:
            s1_char = s1[i - 1]

        if nmmatrix[i][j][DIR] == DIAG:
            ts = s1_char + ts
            bs = s2_char + bs
            if nmmatrix[i][j][MATCH]:
                ms = '|' + ms
            else:
                ms = ' ' + ms

            i = i - 1
            j = j - 1
            continue

        if nmmatrix[i][j][DIR] == LEFT:
            ts = '-' + ts
            bs = s2_char + bs
            ms = ' ' + ms

            j = j - 1
            continue

        if nmmatrix[i][j][DIR] == UP:
            ts = s1_char + ts
            bs = '-' + bs
            ms = ' ' + ms

            i = i - 1
            continue

    return ts, ms, bs


def maxscore(a, b, c):
    A = False
    B = False
    C = False

    if a >= b and a >= c:
        A = True
    if b >= a and b >= c:
        B = True
    if c >= a and c >= b:
        C = True

    return A, B, C


def cost(i, j):
    global gap

    # print(i-1, j-1)
    s1_char = s1[i-1]
    s2_char = s2[j-1]

    # check if diag has been found
    # print(nmmatrix[i-1][j-1], nmmatrix[i-1][j], nmmatrix[i][j-1])

    if nmmatrix[i-1][j-1][SCORE] is None:
        cost(i-1, j-1)

    # check if up has been found
    if nmmatrix[i-1][j][SCORE] is None:
        cost(i - 1, j)

    # check if left has been found
    if nmmatrix[i][j-1][SCORE] is None:
        cost(i, j-1)

    smax = {SCORE: None, DIR: LEFT, MATCH: None}

    leftv = nmmatrix[i][j - 1][SCORE] + gap
    score, match = match_chars(s1_char, s2_char)
    diag = nmmatrix[i - 1][j - 1][SCORE] + score
    upv = nmmatrix[i - 1][j][SCORE] + gap

    A, B, C = maxscore(upv, diag, leftv)

    # take highroad first
    if A:
        smax[SCORE] = upv
        smax[DIR] = UP
    elif B:
        smax[SCORE] = diag
        smax[DIR] = DIAG
        smax[MATCH] = match
    elif C:
        smax[SCORE] = leftv
        smax[DIR] = LEFT

    nmmatrix[i][j][SCORE] = smax[SCORE]
    nmmatrix[i][j][DIR] = smax[DIR]
    nmmatrix[i][j][MATCH] = smax[MATCH]

    # leftv = nmmatrix[i][j - 1][SCORE] + gap
    # smax = {SCORE: leftv, DIR: LEFT, MATCH: None}
    #
    # score, match = match_chars(s1_char, s2_char)
    # diag = nmmatrix[i - 1][j - 1][SCORE] + score
    # if diag > smax[SCORE]:
    #     smax[SCORE] = diag
    #     smax[DIR] = DIAG
    #     smax[MATCH] = match
    #
    # # Find the cost of putting a gap in the second sequence
    # # do this last so that the traceback will be highroad
    # upv = nmmatrix[i-1][j][SCORE] + gap
    # if upv > smax[SCORE]:
    #     smax[SCORE] = upv
    #     smax[DIR] = UP
    #     smax[MATCH] = None
    #
    # nmmatrix[i][j][SCORE] = smax[SCORE]
    # nmmatrix[i][j][DIR] = smax[DIR]
    # nmmatrix[i][j][MATCH] = smax[MATCH]


def alignment_score(ts, bs):
    global gap

    score = 0
    for i in range(0, len(ts)):
        if ts[i] == '-' or bs[i] == '-':
            score = score + gap
        else:
            match_score, _ = match_chars(ts[i], bs[i])
            score = score + match_score

    return score


def do_align(s, m, g, v):
    global s1, s2, gap, nmmatrix, submatrix
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

        s1 = sfile.readline().rstrip().lower()  # read in the first sequence on line 1
        s1_n = len(s1)
        sfile.readline()  # read in the blank line on line 2
        s2 = sfile.readline().rstrip().lower()  # read in the second sequence on line 3
        s2_n = len(s2)

        # Print the first 50 chars of both sequences
        vprint("|\t\tSequence 1 (" + str(s1_n) + "): " + s1[: 50] + "...", v)
        vprint("|\t\tSequence 2 (" + str(s2_n) + "): " + s2[: 50] + "...", v)

    else:
        print("Failed to open " + s, v)
        return

    vprint("|\tDetermining sequence type (DNA/Protein)...", v)

    if 'm' in s1 or 'M' in s1:
        is_protein_seq = True

    vprint("|\t\t-> Protein Sequence" if is_protein_seq else "|\t\t-> DNA Sequence", v)

    if is_protein_seq:
        print("Protein Sequence")
    else:
        print("DNA Sequence")

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
            tup = {'score': None, 'dir': None, 'match': None}
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

    vprint('|\t\tBottom right:', v)
    vprint('|\t\t\tScore ... ' + str(nmmatrix[s1_n][s2_n][SCORE]), v)
    vprint('|\t\t\tDirection ... ' + str(nmmatrix[s1_n][s2_n][DIR]), v)

    vprint('|\n|\tStarting traceback from nmmatrix[' + str(s1_n) + '][' + str(s2_n) + ']', v)
    ts, ms, bs = traceback(s1_n, s2_n)

    print(ts)
    print(ms)
    print(bs)

    print(str(nmmatrix[s1_n][s2_n][SCORE]))



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
    g = args.gpenalty
    verbose = args.verbose
    do_align(sequences, matrix, g, verbose)
