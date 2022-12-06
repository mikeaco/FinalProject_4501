import os
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

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
    return (seq, date, (start_loc, end_loc))

full_list = []
full_list.append("DENV-1_IND_826883_1982.gb")

list_of_files2 = glob.glob("output2"+os.path.sep+"*.gb")
# list_of_files3 = glob.glob("output3"+os.path.sep+"*.gb")
full_list += list_of_files2
writer = open("new_all.fasta", 'w')
for file in full_list:
    seq, date, loc = seq_amino_loc(file)
    writer.write(">"+file + " " + date + "\n"+str(seq)+"\n")

