# need to pip install
# BIO
# blosum
# networkx
    # python -m pip install -U matplotlib
import numpy as np
import blosum as bl
import networkx as nx
import matplotlib.pyplot as plt
from itertools import chain
from collections import defaultdict
import random
from Bio import GenBank
from Bio.Seq import Seq
import sys
import random


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
    j = 1
    for i, char in enumerate(s):
        if char != '-':
            m[i+1] = j
            j+=1
        else:
            m[i+1] = '-'

    return m


s1 = "AAAGAATTCA"
s2 = "AAATGA"
r = showAlignmentG(s1, s2, affineGap, simpleMatch)

data = np.genfromtxt(fname="nextstrain_dengue_denv1_diversity.tsv", delimiter="\t", skip_header=1, filling_values=1)  # change filling_values as req'd to fill in missing values

weights_arr = {}
for i,d in enumerate(data):
    weights_arr[d[0]] = d[1]
print(weights_arr)
counter = {'A':
    {
        'A':112,
        'T':34,
        'G':54,
        'C':12
    },
    'T':
        {
            'A':23,
            'T':113,
            'G':44,
            'C':66
        },
    'G':
        {
            'A':14,
            'T':22,
            'G':142,
            'C':65
        },
    'C':
        {
            'A':33,
            'T':48,
            'G':12,
            'C':111
        }
}

probs = defaultdict(list)

for n in counter:
    total = sum(counter[n].values())
    for m in counter[n]:
        probs[n].append(counter[n][m]/total*100)
print(probs)
def mutate(seq, weight_arr, probs):
    seq = seq.upper()
    chosen_loc = random.choices(
        range(len(seq)), weights=weight_arr, k=1)

    nucl = seq[chosen_loc[0]]
    print(nucl)
    list_prob = probs[nucl]

    nuc_list = ['A','T','G','C']

    new_nucl = random.choices(nuc_list, weights=list_prob, k=1)

    return new_nucl


# testing functionality of mutute
w = [111,23,45,66]
q = 'agct'

print(mutate(q,w,probs))
# random.choice(weights_arr.keys(), p)
years = 1
totmuts = int(.000621 * 1 * list(weights_arr.keys())[-1])

choices = random.choices(list(weights_arr.keys()), list(weights_arr.values()), k =totmuts)
print(choices)
print(weights_arr[choices[0]])
start_loc, end_loc = 0,100



# the_file = sys.argv[1]
from Bio import SeqIO
the_file = "random_sequence.gb"
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

seq, amino, loc = seq_amino_loc(the_file)
seq2, amino2, loc2 = seq_amino_loc("DENV-1_IND_826883_1982.gb")
print(seq)
print(amino)
print(loc)
# map = r[2]
# for i in range(years):
#     choices = random.choices(list(weights_arr.keys()), list(weights_arr.values()), k =totmuts)
#     # print(choices)
#     # print(weights_arr[choices[0]])
#     score = 0
#     for choice in choices:
#         if choice >= loc[0] and choice <= loc[1]:
#             nt = random.choice(['A', 'C', 'G', 'T']) # NEED TO CHANGE THIS TO MIKE
#             seq2[map[choice-start_loc]] = nt


# print(loc, loc2, amino, amino2, seq, seq2, sep="\n")
# r = showAlignmentG(seq, seq2, affineGap, simpleMatch)



