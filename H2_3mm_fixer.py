import csv
from Bio import SeqIO

seq_len = {}
motif_len = {}
motif_type = {}

for record in SeqIO.parse("3mm_seq.fasta", "fasta"):
    seq_len[record.id] = len(record)

with open('nb-lrr_motifs.csv') as data:
    for row in csv.DictReader(data):
        motif_len[row['motif']] = int(row['length'])*3
        motif_type[row['motif']] = row['annotation']

infile = open('H2_3mm.txt', newline='')
csv_in = csv.reader(infile, delimiter='\t',)

with open('H2_3mm_fix.gff', 'w+', newline ='') as outfile:
    csv_out = csv.writer(outfile, delimiter='\t')
    for row in csv_in:
        new_row = row.copy()
        if int(row[4]) == -1:
            row[4] = int(row[3]) + motif_len[row[8][5:]]
        row[2] = motif_type[row[8][5:]]
        csv_out.writerow(row)