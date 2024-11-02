import os

MIN_LENGTH_GEN = 800
DIRECTORY_PATH = "/data/sequeventor"

for filename in os.listdir(DIRECTORY_PATH):
    # Get the full path of the file
    file_path = os.path.join(DIRECTORY_PATH, filename)

    # Check if it's a file (and not a directory)
    if os.path.isfile(file_path):
        print("File:", filename)
    else:
        print("Directory:", DIRECTORY_PATH)
