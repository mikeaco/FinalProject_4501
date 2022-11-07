import time

import numpy as np
from itertools import chain
from Bio import SeqIO
the_file = "random_sequence.gb"
tic = time.perf_counter()
def seq_amino_loc(name):
    amino_acid =""
    seq = ""
    start_loc = 0
    end_loc = 0
    with open(name) as handle:
        record = SeqIO.read(handle, "genbank")
        for i in record.features:
            if i.type == "CDS":
                start_loc = i.location.start
                end_loc = i.location.end
                amino_acid = i.qualifiers["translation"]
        seq = record.seq[start_loc:end_loc]
    return (seq, amino_acid[0], (start_loc, end_loc))


s1, amino, loc = seq_amino_loc(the_file)
s2, amino2, loc2 = seq_amino_loc("DENV-1_IND_826883_1982.gb")
# print(loc, loc2, amino, amino2, seq, seq2, sep="\n")

def alignmentScoreDPG(s1, s2, gapPenalty, match):
    m = np.zeros((len(s1) + 1, len(s2) + 1))
    m[0, 0] = 0
    for i in range(1, len(s1) + 1):
        m[i, 0] = gapPenalty(i)
    for j in range(1, len(s2) + 1):
        m[0, j] = gapPenalty(j)
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            m[i, j] = max(chain((gapPenalty(g) + m[i, j - g] for g in range(1, j + 1)),
                                (gapPenalty(g) + m[i - g, j] for g in range(1, i + 1)),
                                [(match(s1[i - 1], s2[j - 1]) + m[i - 1, j - 1])]))
    return m

def readAlignmentG(s1, s2, m, gapPenalty, match):
    i = len(s1)
    j = len(s2)
    s1a = ""
    s2a = ""
    score = 0
    while i > 0 or j > 0:
        if i > 0 and j > 0 and m[i, j] == m[i - 1, j - 1] + match(s1[i - 1], s2[j - 1]):
            i = i - 1
            j = j - 1
            s1a = s1[i] + s1a
            s2a = (s2[j] if s1[i] == s2[j] else s2[j].upper()) + s2a
            score += match(s1[i], s2[j])
        else:
            foundit = False
            for g in range(1, i + 1):
                #                 if 1 == 1:
                if m[i, j] == m[i - g, j] + gapPenalty(g):
                    s1a = s1[i - g:i] + s1a
                    s2a = ('-' * g) + s2a
                    i = i - g
                    score += gapPenalty(g)
                    foundit = True
                    break
            if not foundit:
                for g in range(1, j + 1):
                    if m[i, j] == m[i, j - g] + gapPenalty(g):
                        s1a = ('-' * g) + s1a
                        s2a = s2[j - g:j] + s2a
                        j = j - g
                        score += gapPenalty(g)
                        foundit = True
                        break
            assert foundit
    return (s1a, s2a, score)

def showAlignmentG(s1, s2, gapPenalty, match):
    m = alignmentScoreDPG(s1, s2, gapPenalty, match)
    #     print(m)
    #     m = /"helo/"
    r = readAlignmentG(s1, s2, m, gapPenalty, match)
    site_map = create_map(r[1])
    print (r[0] + "\n" + r[1] + "\n" + str(r[2]))
    print(site_map)
    return (m, r, site_map)

def affineGap(n, gp = -1, gn = -0.2):
    return gp + (n - 1) * gn

def simpleMatch(a, b):
    return 1 if a == b else -1

def create_map(s):
    m = {}
    j = 0
    for i, char in enumerate(s):
        if char != '-':
            m[i] = j
            j+=1
        else:
            m[i] = '-'

    return m

# s1, a1, loc =
# s1 = "AAAGAAT
# TCA"
# s2 = "AAATGA"
# r = showAlignmentG(seq, seq2, affineGap, simpleMatch)
from Bio import pairwise2
# print('here')
alignments = pairwise2.align.globalxx(s1, s2)
# print(alignments)
# r = showAlignmentG(s1, s2, affineGap, simpleMatch)
# map = {}
counter = {'A':
        {
            'A':0,
            'T':0,
            'G':0,
            'C':0,
            '-':0
        },
    'T':
        {
            'A':0,
            'T':0,
            'G':0,
            'C':0,
            '-':0
        },
    'G':
        {
            'A':0,
            'T':0,
            'G':0,
            'C':0,
            '-':0
        },
    'C':
        {
            'A':0,
            'T':0,
            'G':0,
            'C':0,
            '-':0
        },
    "-":
        {
            'A':0,
            'T':0,
            'G':0,
            'C':0,
            '-':0
        }
}
def create_map2(s1, s2):
    m = {}
    j = 1
    for i, char in enumerate(s1):
        if char != '-' and s2[i] != '-':
            m[i+1] = j
            j+=1
        else:
            m[i+1] = '-'

    return m
ma = create_map2(s1,s2)
# for key,val in ma.items():
#     # print(key, val)
#     if val != "-" and key != "-":
#         # print(key, val)
#         counter[s1[key]][s2[val]] += 1
alignments = sorted(alignments, key= lambda x : x[4])
print(alignments[0][0],alignments[0][1], sep='\n')
a1,a2 = alignments[0][0],alignments[0][1]

for i in range(len(a1)):
    counter[a1[i]][a2[i]] += 1
        # print("ADD"

fil = open("random.txt", 'w')
fil.write(str(counter))
fil.close()
toc = time.perf_counter()
print(f"Time: {toc-tic} seconds")