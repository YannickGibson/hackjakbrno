import subprocess

def bowtie2_align(reference_index, reads_list, path_to_bowtie, output_sam):
    """
    Perform read alignment using Bowtie2 and return the list of starting indices.

    Parameters:
        reference_index (str): Path to the Bowtie2 index files (without file extensions).
        reads_list (list): List of short reads as strings to align.
        output_sam (str): Path to the output SAM file to save the alignment results.

    Returns:
        list: A list containing the starting indices of each read in the alignment.
    """

    # Create a temporary interleaved FASTQ file to store the reads
    temp_fastq = "temp_interleaved.fastq"
    with open(temp_fastq, "w") as f:
        for i in range(0, len(reads_list)):
            f.write(f"@read_{i}\n{reads_list[i]}\n+\n{'I'*len(reads_list[i])}\n")  # Assuming all reads have the same length

    # Command to run Bowtie2 alignment with interleaved input
    bowtie2_command = f"{path_to_bowtie}/bowtie2 -x {reference_index} -q {temp_fastq} -S {output_sam}" # --interleaved"

    try:
        # Execute Bowtie2 command using subprocess
        subprocess.run(bowtie2_command, shell=True, check=True)
        print("Bowtie2 alignment completed successfully.")
        # Parse the SAM file to extract the starting indices of each read
        starting_indices = []
        with open(output_sam, "r") as sam_file:
            for line in sam_file:
                if not line.startswith("@"):  # Skip header lines
                    fields = line.strip().split("\t")
                    if len(fields) >= 4:
                        starting_index = int(fields[3])
                        starting_indices.append(starting_index)

        print("Bowtie2 alignment completed successfully.")
        return starting_indices

    except subprocess.CalledProcessError as e:
        print(f"Error during Bowtie2 alignment: {e}")

    # Remove the temporary interleaved FASTQ file
    subprocess.run(f"rm {temp_fastq}", shell=True, check=True)
    
    starting_indices = []
    with open(output_sam, "r") as sam_file:
        for line in sam_file:
            if not line.startswith("@"):  # Skip header lines
                fields = line.strip().split("\t")
                if len(fields) >= 4:
                    starting_index = int(fields[3])
                    starting_indices.append(starting_index)
    
    return starting_indices