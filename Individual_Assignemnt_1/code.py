#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 19:23:42 2024

@author: inestavares
"""

def create_submat (match, mismatch, alphabet):
    """ substitution matrix as dictionary """
    sm = {}
    for c1 in alphabet:
        for c2 in alphabet:
            if (c1 == c2):
                sm[c1+c2] = match
            else:
                sm[c1+c2] = mismatch
    return sm

# score of a position (column)
def score_pos(c1, c2, sm, g):
    """score of a position (column)"""
    if c1 == "-" or c2=="-":
        return g
    else:
        return sm[c1+c2]


def score_affinegap (seq1, seq2, sm, g, r):
    ''' calculates the score of alignment based on : affine_gap(len) = g + r*len
    if the gap is open (first occurrence) sum value g; if gap continues sum r to each new gap position;
    if there is no gap use the substitution matrix for the score.
    '''
    res = 0
    ingap1 = False # two f are true when inside gap sequences
    ingap2 = False
    for i in range(len(seq1)):
        if seq1[i]=="-":
            # gap is already open; add r
            if ingap1: res += r
            else:
                # gap is open for the first time; add g
                ingap1 = True
                res += g
        elif seq2[i]=="-":
            # gap is already open; add r
            if ingap2: res += r
            else:
                # gap is open for the first time; add g
                ingap2 = True
                res += g 
        else:
            # no gaps; use substitution matrix
            if ingap1: ingap1 = False
            if ingap2: ingap2 = False
            res += sm[seq1[i]+seq2[i]]
    return res


## global alignment 

def needleman_Wunsch (seq1, seq2, sm, g):
    """Global Alignment"""
    S = [[0]]
    T = [[0]]
    # initialize gaps in rows
    for j in range(1, len(seq2)+1):
        S[0].append(g * j)
        T[0].append(3)  # horizontal move: 3
    # initialize gaps in cols
    for i in range(1, len(seq1)+1):
        S.append([g * i])
        T.append([2])  # vertical move: 2
    # apply the recurrence to fill the matrices
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos (seq1[i], seq2[j], sm, g); # diagonal
            s2 = S[i][j+1] + g  # vertical
            s3 = S[i+1][j] + g # horizontal
            S[i+1].append(max(s1, s2, s3)) # na matrix score add max value
            T[i+1].append(max3t(s1, s2, s3))
    return (S, T)

def max3t (v1, v2, v3):
    """Provides the integer to fill in T"""
    if v1 > v2:
        if v1 > v3: return 1
        else: return 3
    else:
        if v2 > v3: return 2
        else: return 3

def recover_align (T, seq1, seq2):
    # alignment are two strings
    res = ["", ""]
    i = len(seq1)
    j = len(seq2)
    while i>0 or j>0:
        if T[i][j]==1: # diagonal move
            res[0] = seq1[i-1] + res[0] # add to align of seq1 a symbol from seq1(i-1)
            res[1] = seq2[j-1] + res[1] # add to align of seq2 a symbol from seq2(i-1)
            i -= 1
            j -= 1
        elif T[i][j] == 3: # horizontal move
            res[0] = "-" + res[0]   # insert gap na seq 1
            res[1] = seq2[j-1] + res[1] # insert symbol from seq2
            j -= 1
        else: # vertical move
            res[0] = seq1[i-1] + res[0] # insert symbol from seq1
            res[1] = "-" + res[1] # insert gap na seq 2
            i -= 1
    return res

## local alignment

def smith_Waterman (seq1, seq2, sm, g):
    """Local alignment"""
    S = [[0]]
    T = [[0]]
    maxscore = 0
    # first row filled with zero
    for j in range(1, len(seq2)+1):
        S[0].append(0)
        T[0].append(0)
    # first column filled with zero
    for i in range(1, len(seq1)+1):
        S.append([0])
        T.append([0])
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos (seq1[i], seq2[j], sm, g); 
            s2 = S[i][j+1] + g
            s3 = S[i+1][j] + g
            b = max(s1, s2, s3)
            if b <= 0:
                S[i+1].append(0)
                T[i+1].append(0)
            else:
                S[i+1].append(b)
                T[i+1].append(max3t(s1, s2, s3))
                if b > maxscore: 
                    maxscore = b
    return (S, T, maxscore)

def recover_align_local (S, T, seq1, seq2):
    """recover one of the optimal alignments"""
    res = ["", ""]
    """determine the cell with max score"""
    i, j = max_mat(S)
    """terminates when finds a cell with zero"""
    while T[i][j]>0:
        if T[i][j]==1:
            res[0] = seq1[i-1] + res[0]
            res[1] = seq2[j-1] + res[1]
            i -= 1
            j -= 1
        elif T[i][j] == 3:
            res[0] = "-" + res[0];
            res[1] = seq2[j-1] + res[1] 
            j -= 1
        elif T[i][j] == 2:
            res[0] = seq1[i-1] + res[0]
            res[1] = "-" + res[1]
            i -= 1
    return res

def coords_appendix(coords_mat, row, col):
    """Append coordinates to the list and return the updated list"""
    resmat = (row, col)
    if coords_mat is None:
        coords_mat = []  # Initialize the list if it's None
    coords_mat.append(resmat)
    return coords_mat

def max_mat(mat):
    """finds the max cell in the matrix"""
    maxval = mat[0][0]
    maxcoords = [(0, 0)]
    nrows = len(mat)
    ncols = len(mat[0])
    # returns the cell with maximum value
    for r in range(nrows):
        for c in range(ncols):
            if mat[r][c] > maxval:
                maxval = mat[r][c]
                maxcoords = []
                maxcoords = coords_appendix(maxcoords, r, c)
            elif mat[r][c] == maxval:
                maxcoords = coords_appendix(maxcoords, r, c)
            else:
                continue
    return maxcoords[0]


def print_mat (mat):
    for i in range(0, len(mat)):
        print(mat[i]) 

def global_and_local_alignment(s1,s2,match=3,mismatch=-1,gap=-3):
    #Function to obtain the global and local alignment
    
    #Global alignment
    global_st = needleman_Wunsch(s1, s2, 
                                 create_submat(match, mismatch, "ACGT"), gap)
    global_S = global_st[0]  # Score matrix
    global_T = global_st[1]  # Traceback matrix
    global_bs=global_bs = global_S[len(s1)][len(s2)]
    global_alignment = recover_align(global_T, s1, s2)  # Optimal alignment


    # Local Alignment
    local_st = smith_Waterman(s1, s2, 
                              create_submat(match, mismatch, "ACGT"), gap)
    local_S = local_st[0]  # Score matrix
    local_T = local_st[1]  # Traceback matrix
    local_bs = local_st[2]  # Best score
    local_alignment = recover_align_local(local_S, local_T, s1, s2)  # Optimal alignment
    return ("Global", global_S, global_T, global_bs, global_alignment),("Local", local_S, local_T, local_bs, local_alignment)


# Test the alignment_info function

s1 = "CGCTTA"
s2 = "ACCTAA"
match = 3
mismatch = -1
gap = -3
    
global_info,local_info = global_and_local_alignment(s1, s2, match, mismatch, gap)

def print_alignment_info(info):
    align_type, sm, TM, bs, opt_align = info
    print("Alignment Type:", align_type)
    print("Score Matrix:")
    for i in sm:
        print(i)
    print("Traceback Matrix:")
    for i in TM:
        print(i)
    print("Best Score:", bs)
    print("Optimal Alignment:")
    for line in opt_align:
        print(line)

print("Global Alignment:")
print_alignment_info(global_info)

print("\nLocal Alignment:")
print_alignment_info(local_info)
