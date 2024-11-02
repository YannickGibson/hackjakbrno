from model import DNABERT2
from Bio import SeqIO
from tqdm import tqdm
import os
import gzip
import pandas as pd

SEQUENCING_DATA_PATH = "data/sequencing_data".replace('\\', "/")
EXPORT_EVERY = 10
EXPORT_PATH = "data/rows_for_db.csv"


def load_fastq_file(fastq_path: str, model: DNABERT2, rows: list) -> None:
    
    with gzip.open(fastq_path, 'rt') as file:
        for i, record in tqdm(enumerate(SeqIO.parse(file, 'fastq')), total=4000):
            # Convert sequence to embedding
            vector = model.get_embedding(str(record.seq)).detach().tolist()[0]

            rows.append([fastq_path, i, vector])
            if i % EXPORT_EVERY == 0:
                print(f"Exporting {len(rows)} rows to '{EXPORT_PATH}'.")
                df = pd.DataFrame(rows, columns=["file_path", "index", "vector"])
                df.to_csv(EXPORT_PATH, index=False, sep=";")


def load_barcode(folder_path: str, model: DNABERT2, rows: list) -> None:
    fastq_files = sorted(os.listdir(folder_path))
    fastq_files = [file for file in fastq_files if file.endswith(".fastq.gz")]
    print(f"Opening {len(fastq_files)} fastq files.")
    for fastq_file in fastq_files:
        load_fastq_file(fastq_path=folder_path + "/" + fastq_file, model=model, rows=rows)


def main() -> None:
    # Load model
    model = DNABERT2()
    #model = None

    # create df to append rows to 
    rows = []  # file_path, index, vector

    # Load barcode folders
    barcode_folders = sorted(os.listdir(SEQUENCING_DATA_PATH))
    print(f"Exporting {len(barcode_folders)} barcode folders.")
    for folder_name in barcode_folders:
        load_barcode(folder_path=SEQUENCING_DATA_PATH + "/" + folder_name, model=model, rows=rows)
        

if __name__ == '__main__':
    main()
