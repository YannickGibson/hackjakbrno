from model import DNABERT2
from Bio import SeqIO
from tqdm import tqdm
import os
import gzip
import pandas as pd
from natsort import natsorted

SEQUENCING_DATA_PATH = "data/sequencing_data".replace('\\', "/")
EXPORT_EVERY = 1000
EXPORT_PATH = "data/rows_for_db.csv"
MAX_SEQUENCE_LENGTH = 10000


def load_fastq_file(fastq_path: str, model: DNABERT2, rows: list) -> None:
    
    with gzip.open(fastq_path, 'rt') as file:
        pbar = tqdm(enumerate(SeqIO.parse(file, 'fastq')), total=4000)
        for i, record in pbar:
            # Convert sequence to embedding
            sequence_string = str(record.seq)

            # Truncate string if too long
            sequence_string = sequence_string[:MAX_SEQUENCE_LENGTH]

            # Update progress bar info before inference
            pbar.set_description(f"Sequence length: {len(sequence_string)}, fastq_path: '{fastq_path}'")

            # Inference: get embedding
            vector = model.get_embedding(sequence_string).detach().tolist()[0]

            rows.append([fastq_path, i, sequence_string, ",".join(str(num) for num in vector)])
            # show sequence_string len
            if i == 0 or (i - 1) % EXPORT_EVERY == 0 or i == 3999:  # first or every x'th or last
                print(f"Exporting {len(rows)} rows to '{EXPORT_PATH}'.")
                df = pd.DataFrame(rows, columns=["file_path", "index", "sequence_string", "sequence_vector"])
                df.to_csv(EXPORT_PATH, index=False, sep=";")


def load_barcode(folder_path: str, model: DNABERT2, rows: list) -> None:
    fastq_files = natsorted(os.listdir(folder_path))
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
