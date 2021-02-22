import csv
from Bio import SeqIO

# Protip: Don't initialise massive contigs as a dictionary as it will kill windows... Instead iterate through each time (less memory)

def sortString(str):
    str = str.replace(',', '')
    return ''.join(sorted(str)) 

def StringIndexReplace (str, index, replacement):
    new_str = str
    if index < len(str):
        new_str = str[0:index] + replacement + str[index + 1:]
    return new_str

# Dictionary to store parent SNPs as "unitigname_SNPposition: IUPACcode"
parent_SNPs = {}

# List of the unique unitigs present in commonSNPs, for sorting parental SNPs
unique_unitigs = []

# IPUAC dictionary for determining SNP type
IUPAC = {
    'AG' : 'R',
    'CT' : 'Y',
    'CG' : 'S',
    'AT' : 'W',
    'GT' : 'K',
    'AC' : 'M',
    'CGT' : 'B',
    'AGT' : 'D',
    'ACT' : 'H',
    'ACG' : 'V',
}

# Get the list of unitigs we care about, for later filtering
with open('commonSNPs5.sort.fixed.vcf') as infile:

        for row in csv.DictReader(infile, delimiter='\t'):

            if row['#CHROM'] not in unique_unitigs:

                unique_unitigs.append(row['#CHROM'])

# Clean up parental SNPs so only unitg SNPs we care about are carried over
with open('parent_SNPs_popped.vcf', 'w', newline ='') as outfile, open('parent_SNPs.vcf') as infile:

    for line in infile:

        # JANKY, DON'T TRUST THIS METHOD FOR EXPANSION.
        # ACTUALLY, DON'T TRUST ANY OF THIS.
        if  line[:6] == '#CHROM' or line[:10] in unique_unitigs:

            outfile.write(line)

# Generate the dictionary of parental SNPs, nasty string manipulation
with open('parent_SNPs_popped.vcf') as infile:

    for row in csv.DictReader(infile, delimiter='\t'):

        row_code = row['REF'] + row['ALT']
        row_code = sortString(row_code)

        parent_SNPs[row['#CHROM'] + '_' + str(row['POS'])] = IUPAC[row_code]

# Generate formatted strings and output as tsv
with open('commonSNPs5.sort.fixed.vcf') as infile, open('KASP_markers.csv', 'w', newline='') as outfile:

    csv_writer = csv.writer(outfile)
    csv_writer.writerow(['KASP_Marker', 'KASP_Sequence', 'SNPs'])

    for row in csv.DictReader(infile, delimiter='\t'):

        for record in SeqIO.parse("3mm_complete.fasta", "fasta"):

            if record.id == row['#CHROM']:

                KASP_Marker = record.id + '_' + row['POS']

                basepos = int(row['POS']) - 1 # -1 as index from 0

                KASP_Sequence = str(record.seq[basepos-50:basepos+51])

                SNP_Count = 0

                for i in range(basepos-50, basepos+51): # A very likely error will occur when an SNP is incident within 50bp of a unitig end

                    native = i - basepos + 50
                    
                    key = record.id + '_' + str(i+1) # This makes me laugh every time

                    if key in parent_SNPs:

                        KASP_Sequence = StringIndexReplace(KASP_Sequence, native, parent_SNPs[key])

                        SNP_Count +=1

                KASP_Sequence = KASP_Sequence[:50] + '[' + KASP_Sequence[50] + ']' + KASP_Sequence[51:]

                csv_writer.writerow([KASP_Marker, KASP_Sequence, SNP_Count])