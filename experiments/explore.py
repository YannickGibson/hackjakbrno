import pandas as pd

LABELS_PATH = 'data/KP_samples_with_resistance_tab.xlsx'

# convert excel to pandas dataframe
df = pd.read_excel(LABELS_PATH, header=1)

df.head()

# descend to parent dir
import sys
import os
__file__ = os.path.abspath("")
sys.path.insert(1, os.path.abspath( os.path.dirname(__file__)) )

import pandas as pd
from solver import solve, mock_solve

# Choose mock method (optional)
solve = solve

n_similar_vectors = [2000]
avg_accuracies = []  # avg per every barcode for each hyper parameter
results = []

def get_accuracy(result: dict[dict[str: bool]], df: pd.DataFrame) -> float:
    accuracies = []
    for barcode_name, bacterias in result.items():
        correct_fields = []
        for bacteria_name, detected in bacterias.items():
            label = df[df['barcode'] == barcode_name][bacteria_name].values[0]
            is_correct = label == detected
            correct_fields.append(is_correct)
        accuracies.append(sum(correct_fields) / len(correct_fields))
    return accuracies

def get_avg_accuracy(result: dict[dict[str: bool]], df: pd.DataFrame) -> float:
    accuracies = get_accuracy(result, df)
    return sum(accuracies) / len(accuracies)

for n in n_similar_vectors:
    # Solve
    result = solve(n, verbose=False)
    results.append(result)
    # Calculate avg accuracy for each barcode
    avg_accuracy = get_avg_accuracy(result, df)
    avg_accuracies.append(avg_accuracy)

print(f"Average accuracies: {avg_accuracies}")

import matplotlib.pyplot as plt

# plot line graph
plt.plot(n_similar_vectors, avg_accuracies, marker='x')

# set labels
plt.xlabel('Number of best similaries')
plt.ylabel('Mean Accuracy')

plt.show()