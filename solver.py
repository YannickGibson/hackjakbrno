import os
from Bio import SeqIO
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd
from natsort import natsorted

from utils.check import is_match_pairwise2, is_match_waterman
from model.model import DNABERT2
from modelRDE import RDE


def solve(results_number: int, rcomplement: bool = False, plot_gene_scores: bool = False, verbose=True) -> dict[str, dict[str, bool]]:
    from iris_database import IrisDatabaseHandler

    # Arguments
    barcodes_path = "data/sequencing_data".replace("\\", "/")
    genes_path = "data/resistance_genes_sequence".replace("\\", "/")

    # model = RDE()
    model = DNABERT2()
    model_str = "DNABERT2"

    # Load genes
    genes = []
    genes_vectors = {}
    gene_names = []
    gene_files = natsorted(os.listdir(genes_path))
    for i, file_name in enumerate(tqdm(gene_files)):
        with open(genes_path + "/" + file_name, 'r') as gene_file:
            gene_names.append(".".join(file_name.split(".")[:-1]))
            genes.append([])
            genes_vectors[i] = []
            # Append all records: gene and rcomplement
            for record in SeqIO.parse(gene_file, 'fasta'):
                genes[-1].append(record)
                # genes_vectors[i].append(model.get_embedding(str(record.seq)).tolist())
                genes_vectors[i].append(model.get_embedding(str(record.seq)).detach().tolist()[0])
                if not rcomplement:
                    break

    # Init database
    database = IrisDatabaseHandler(results_number=results_number)

    # Load barcodes
    barcodes: list[str] = [item[0] for item in database.get_unique_barcodes()]
    print(barcodes)# ["barcode20"]  # os.listdir(barcodes_path)
    print(database.get_unique_barcodes())
    print(f"Testing barcodes: '{barcodes}' for the following genes: '{gene_names}'")
    barcode_gene_positive = {barcode: [False] * len(genes) for barcode in barcodes}

    # Check barcodes
    sequenator_finished = True
    while True:
        iterable = barcodes if verbose else tqdm(barcodes)
        for barcode_name in iterable:
            if verbose:
                print(f"Checking barcode: {barcode_name}")
            # Check if barcode is in genes
            scores = []
            for gene_index, gene in enumerate(tqdm(genes)):
                if verbose:
                    print(f"Current gene: {gene_names[gene_index]}")
                if barcode_gene_positive[barcode_name][gene_index]:  # constraint is already satisfied
                    continue
                gene_list_str = str(genes_vectors[gene_index][0])  # get string sequence

                # Cosine similarity vector database search
                matches: list[tuple] = database.search(barcode_name, gene_list_str, model_str)
                if verbose:
                    print(f"Checking {len(matches)} amount of matches.")

                # Verify similar vectors
                for barcode_name, file_path, index, sequence_string, sequence_vector, model_name in matches:
                    for curr_gene in gene:  # gene or its rcomplement
                        is_match, score = is_match_waterman(str(curr_gene.seq), sequence_string, verbose=verbose, verbose_if_matched=True)
                        scores.append(score)
                        if is_match:
                            if verbose:
                                print(f"Filename of found bacteria: {file_path}\nIndex: {index}")
                            barcode_gene_positive[barcode_name][gene_index] = True
                            break

                # Set up the figure and axis
                if plot_gene_scores:
                    # Init plot
                    fig, ax = plt.subplots()
                    # Plot
                    ax.plot(range(0, len(scores)), scores, lw=2, label="Scores")
                    ax.plot(range(0, len(scores)), pd.Series(scores).rolling(window=15).mean().tolist(), lw=2, label="Rolling Mean")
                    # Configure
                    ax.set_xlabel("Vector Index")
                    ax.set_ylabel("Scores")
                    # Display
                    plt.show(block=True)

        # End if sequenator is done
        if sequenator_finished:
            break
        else:
            pass  # sleep for x seconds

    result = {barcode_name: {}}
    # Print results
    for barcode_name, is_gene_positive_list in barcode_gene_positive.items():
        for gene_index, is_gene_positive in enumerate(is_gene_positive_list):
            gene_name = gene_names[gene_index]
            print(f"SSSS: {gene_name}")
            print(f"SSSS: {barcode_name}")
            result[barcode_name][gene_name] = is_gene_positive
            print(f"Barcode: {barcode_name}, Gene: {gene_name}, Positive: {is_gene_positive}")
        print("\n")

    return result


def mock_solve(results_number: int, rcomplement: bool = False, plot_gene_scores: bool = False, verbose=True) -> dict[str, dict[str, bool]]:
    # Mock results when the database is not available. Higher results number yields better results.
    return {
        "barcode20": {
            "aph(6)-Id": True    if results_number >= 100 else False,  # Resistance 0
            "aph(3'')-Ib": True  if results_number >= 80 else False,  # Resistance 0
            "aac(3)-IIa": True   if results_number >= 50 else False,  # Resistance 1
            "aac(3)-IId": False  if results_number >= 10 else True,  # Resistance 1
            "tet(A)": True       if results_number >= 5 else False,  # Resistance 2
            "tet(D)": False      if results_number >= 2 else True,  # Resistance 2
        }
    }


def main() -> None:
    solve(
        results_number=10,  # Number of results for cosine similarity search db to search through
        rcomplement=False
    )


if __name__ == "__main__":
    main()

