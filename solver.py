import os
from Bio import SeqIO
from tqdm import tqdm
from model import DNABERT2

from iris_database import IrisDatabaseHandler
from check import is_match_pairwise2

def main() -> None:

    # Arguments
    barcodes_path = "data/sequencing_data".replace("\\", "/")
    genes_path = "data/resistance_genes_sequence".replace("\\", "/")

    model = DNABERT2()

    # Load genes
    genes = []
    genes_vectors = {}
    gene_names = []
    files = os.listdir(genes_path)
    for i, file in enumerate(tqdm(files)):
        with open(genes_path + "/" + file, 'r') as gene_file:
            gene_names.append(".".join(file.split(".")[:-1]))
            genes.append([])
            genes_vectors[i] = []
            # Append all records: gene and rcomplement
            for record in SeqIO.parse(gene_file, 'fasta'):
                genes[-1].append(record)
                genes_vectors[i].append(model.get_embedding(str(record.seq)).detach().tolist()[0])

    # Init database
    database = IrisDatabaseHandler()

    # Load barcodes
    barcodes = ["barcode20"]  # os.listdir(barcodes_path)
    barcode_gene_positive = {barcode: [False] * len(genes) for barcode in barcodes}

    # Check barcodes
    sequenator_finished = True
    while True:
        for barcode_name in barcodes:
            print(f"Checking barcode: {barcode_name}")
            # Check if barcode is in genes
            for gene_index, gene in enumerate(tqdm(genes)):
                if barcode_gene_positive[barcode_name][gene_index]:  # constraint is already satisfied
                    continue
                gene_list_str = str(genes_vectors[gene_index][0])  # get string sequence

                # Cosine similarity vector database search
                matches: list[tuple] = database.search(barcode_name, gene_list_str)
                print(f"Checking {len(matches)} amount of matches.")

                # Verify similar vectors
                for barcode_name, file_path, index, sequence_string, sequence_vector in matches:
                    is_match, score = is_match_pairwise2(str(gene[0].seq), sequence_string, verbose=True)
                    # print(f"IS MATCH: {is_match}  AND SCORE: {score}")
                    if is_match:
                        barcode_gene_positive[barcode_name][gene_index] = True
                        break
        
        # End if sequenator is done
        if sequenator_finished:
            break

    # Print results
    for barcode_name, is_gene_positive_list in barcode_gene_positive.items():
        for gene_index, is_gene_positive in enumerate(is_gene_positive_list):
            gene_name = gene_names[gene_index]
            print(f"Barcode: {barcode_name}, Gene: {gene_name}, Positive: {is_gene_positive}")
        print("\n")


if __name__ == "__main__":
    main()
