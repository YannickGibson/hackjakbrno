from model.model import DNABERT2
from Bio import SeqIO
from tqdm import tqdm
import os
import pandas as pd
from iris_database import IrisDatabaseHandler

MIN_LENGTH_GEN=800
MAX_SEQUENCE_LENGTH=10000
FASTQ_DNA_SEQUENCE_PATH="data/sequeventor/"
MODEL=DNABERT2()
DB = IrisDatabaseHandler()

for filename in os.listdir(FASTQ_DNA_SEQUENCE_PATH):
    if not filename.endswith(".fastq"):
        print("Current pointer is not an actual file!")
        continue

    # Get the full path of the file
    file_path = os.path.join(FASTQ_DNA_SEQUENCE_PATH, filename)

    with open(file_path, 'r') as file:
        pbar = tqdm(enumerate(SeqIO.parse(file, 'fastq')), total=4000)
        for i, record in pbar:
            # Convert sequence to embedding
            sequence_string = str(record.seq)
            print(sequence_string)

            # Truncate string if too long
            sequence_string = sequence_string[:MAX_SEQUENCE_LENGTH]

            # Update progress bar info before inference
            pbar.set_description(f"Sequence length: {len(sequence_string)}")

            # Inference: get embedding
            vector = MODEL.get_embedding(sequence_string).detach().tolist()[0]

            rows = [[file_path.split("_")[2], file_path, i, sequence_string, ",".join(str(num) for num in vector)]]
            df = pd.DataFrame(rows, columns=["barcode", "file_path", "index", "sequence_string", "sequence_vector"])

            DB.insert_embedded_sequence(df)
