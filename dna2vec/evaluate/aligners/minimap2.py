import subprocess

def minimap2_align(reference_index, reads_list, path_to_minimap_folder, output_sam):
    """
    Perform read alignment using minimap2 and return the list of starting indices.

    Parameters:
        reference_index (str): Path to the minimap2 index files (without file extensions).
        reads_list (list): List of short reads as strings to align.
        output_sam (str): Path to the output SAM file to save the alignment results.

    Returns:
        list: A list containing the starting indices of each read in the alignment.
    """
    # Create a temporary FASTA file to store the reference genome
    # temp_reference = "temp_reference.fasta"
    # with open(temp_reference, "w") as f:
    #     f.write(">Reference\n" + reference_index + "\n")

    # Create a temporary FASTQ file to store the reads
    temp_fastq = "temp_reads.fastq"
    with open(temp_fastq, "w") as f:
        for i, read in enumerate(reads_list):
            f.write(f"@read_{i}\n{read}\n+\n{'I'*len(read)}\n")  # Assuming all reads have the same length

    # Command to run minimap2 alignment
    minimap2_command = f"{path_to_minimap_folder}/minimap2/minimap2 -ax sr {reference_index} {temp_fastq} > {output_sam}"

    try:
        # Execute minimap2 command using subprocess
        subprocess.run(minimap2_command, shell=True, check=True)
        print("minimap2 alignment completed successfully.")

        # Parse the SAM file to extract the starting indices of each read
        starting_indices = []
        with open(output_sam, "r") as sam_file:
            for line in sam_file:
                if not line.startswith("@"):  # Skip header lines
                    fields = line.strip().split("\t")
                    if len(fields) >= 4:
                        starting_index = int(fields[3])
                        starting_indices.append(starting_index)

        return starting_indices

    except subprocess.CalledProcessError as e:
        print(f"Error during minimap2 alignment: {e}")

    finally:
        # Remove the temporary files
        subprocess.run(f"rm {temp_fastq}", shell=True, check=True)
