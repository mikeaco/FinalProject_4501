from itertools import chain
import numpy as np
import blosum as bl
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
import math
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import defaultdict
import glob
import os

pre_set = {'A': {'A': 0, 'T': 7, 'G': 3, 'C': 6}, 'T': {'A': 8, 'T': 0, 'G': 1, 'C': 8}, 'G': {'A': 17, 'T': 6, 'G': 0, 'C': 5}, 'C': {'A': 3, 'T': 20, 'G': 5, 'C': 0}}
# pre_set = {'A': {'A': 13025, 'T': 7, 'G': 3, 'C': 6}, 'T': {'A': 8, 'T': 8030, 'G': 1, 'C': 8}, 'G': {'A': 17, 'T': 6, 'G': 9955, 'C': 5}, 'C': {'A': 3, 'T': 20, 'G': 5, 'C': 7275}}

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

def probability_from_dict(dict):
    probs = defaultdict(list)

    for n in dict:
        total = sum(dict[n].values())
        for m in dict[n]:
            probs[n].append(dict[n][m]/total*100)
    return probs

def mutate(probs, nucleotide):
    list_prob = probs[nucleotide]
    nuc_list = ['A','T','G','C']
    new_nucl = random.choices(nuc_list, weights=list_prob, k=1)
    return new_nucl

def seq_amino_loc(name):
    amino_acid =""
    # seq = ""
    start_loc = 0
    end_loc = 0
    # date = ""
    with open(name) as handle:
        record = SeqIO.read(handle, "genbank")
        date = record.annotations['date']
        for j in record.features:
            if j.type == "CDS":
                start_loc = j.location.start
                end_loc = j.location.end
                amino_acid = j.qualifiers["translation"]
        seq = record.seq[start_loc:end_loc]
    return (seq, amino_acid[0], (start_loc, end_loc, date))

def create_map(s1, s2):
    m = {}
    j = 0
    for i, char in enumerate(s1):
        if char != '-' and s2[i] != '-':
            m[i] = j
            j+=1
        else:
            x = ''
            if char == '-':
                x = j+1
                j+=2
            if s2[i] == '-':
                x = i
            m[i] = x
    return m



parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('origin_file', type = str, help='genbank file with corrispondence to diversity file', default="DENV-1_IND_826883_1982.gb")
parser.add_argument('diversity_file', type = str, help='diversity file', default="nextstrain_dengue_denv1_diversity.tsv")
parser.add_argument('mutation_rate', type = float, help='mutation rate', default=.000621)
parser.add_argument('mutation_time', type = float, help='number of years to mutation', default=1.0)
parser.add_argument('path', type=str, help='Path to genbank files', default="genbank_files")
parser.add_argument('save_path', type=str, help='Desired save location', default=".")
args = parser.parse_args()
list_of_files = glob.glob(args.path+os.path.sep+"*.gb")
# for i in
print("Reading diversity file...", end='')
data = np.genfromtxt(fname=args.diversity_file, delimiter="\t", skip_header=1, filling_values=1)  # change filling_values as req'd to fill in missing values
weights_arr = {}
for i, d in enumerate(data):
    weights_arr[d[0]] = d[1]
print("Finished")
seq, amino, loc = seq_amino_loc(args.origin_file)
file_num = 0
for gb_file in list_of_files:
    change_file = gb_file
    file_num += 1
    for year in range(math.ceil(args.mutation_time)):
        print("Reading in genbank files and extracting gene sequence and locations...", end='')
        seq2, amino2, loc2 = seq_amino_loc(change_file)
        print("Finished")

        print("Aligning sequences and choosing best alignment (takes < 1 minute for dengue data)...", end='')
        alignments = pairwise2.align.globalxx(seq, seq2)
        alignments = sorted(alignments, key= lambda x : x[4])
        a1, a2 = alignments[0][0], alignments[0][1]
        print("Finished")

        print("Creating map...", end='')
        the_map = create_map(a1, a2) # the map implementation
        print("Finished")

        totmuts = math.ceil(args.mutation_rate * list(weights_arr.keys())[-1])
        blosum_matrix = bl.BLOSUM(62)
        seq2 = str(seq2)
        score = 0.0
        choices = random.choices(list(weights_arr.keys()), list(weights_arr.values()), k =totmuts)

        for choice in choices:
            if choice >= loc[0] and choice <= loc[1]:
                # nt = random.choice(['A', 'C', 'G', 'T']) # NEED TO CHANGE THIS TO MIKES Markov Model
                locally = the_map[choice-loc[0]]
                probs = probability_from_dict(pre_set)
                nt = mutate(probs, seq2[locally])
                amino_before = Seq(seq2[locally - (locally%3): locally - (locally%3) +3]).translate()
                seq2 = seq2[:locally] + nt[0] + seq2[locally+1:]
                amino_after = Seq(seq2[locally - (locally%3): locally - (locally%3) +3]).translate()
                score += blosum_matrix[amino_before+amino_after]
                if amino_before != amino_after:
                    amino_change_loc = locally//3
                    amino2 = str(amino2)
                    amino2 = amino2[:amino_change_loc] + amino_after + amino2[amino_change_loc+1:]
        amino2 = str(amino2)
        seq_str = Seq(seq2)
        record = SeqRecord(seq=seq_str, id=str(math.ceil(args.mutation_time)), name="Changed_file", description="Mutated file with a score of " + str(score) +" "+ str(math.ceil(args.mutation_time))+" year after the original")
        feature = SeqFeature(FeatureLocation(start=0, end=(len(seq2)+1)), type='CDS', qualifiers={'translation': [amino2]})
        changed_date = loc2[2]
        changed_date = int(changed_date[-4:]) + 1
        record.annotations['date'] = loc2[2][:-4] + str(changed_date)
        record.features.append(feature)
        filename = "updated"+str(file_num)+".gb"
        complete_path = os.path.join(args.save_path, filename)
        output_file = open(filename, "w")
        SeqIO.write(record, complete_path, "genbank")
        change_file = complete_path
        output_file.close()
        print(f"COMPLETED {year+1} year(s)")

print("Program ended successfully")
