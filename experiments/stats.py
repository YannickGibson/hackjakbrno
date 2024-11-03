import os
from Bio import SeqIO

#GENE_FOLDER = r'/home/user_pool_2/hackathon/resistance_genes_sequence/'
SEQUENCE_FOLDER = r'/home/ivan/school/hackjakbrno/barcode25/'
# load all files in the GENE_FOLDER

def count_stats(folder, verbose=False, gzip=False, extension="fasta"):
    min_len = 99999999
    avg_len = 0
    max_len = 0
    count = 0

    for file in os.listdir(folder):
        with open(folder + file, 'r') as gene_file:
            print(f"File: {gene_file}")
            for record in SeqIO.parse(gene_file, extension):
                count += 1
                l = len(record)
                if verbose:
                    print(f"[{count}] lenght: {l}")
                avg_len += l
                if min_len > l:
                    min_len = l
                if max_len < l:
                    max_len = l
    avg_len /= count
    print(f"Folder: {folder}\nAvg: {avg_len}, min: {min_len}, max: {max_len}, count: {count}")

#count_stats(GENE_FOLDER, verbose=True, extension="fasta")
count_stats(SEQUENCE_FOLDER, extension="fastq")