import re
import os
import statistics

def parse_log_file(file_path):
    nni_data = []
    optimal_score = None
    total_runtime = 0
    total_trees = 0

    with open(file_path, 'r') as file:
        current_nni = {}
        for line in file:
            if line.startswith("NNI Iter:"):
                if current_nni:
                    nni_data.append(current_nni)
                current_nni = {}
            elif "Number of trees checked:" in line:
                current_nni['trees'] = int(line.split(":")[1].strip())
            elif "Current score:" in line:
                current_nni['score'] = float(line.split(":")[1].strip())
            elif "Runtime (s):" in line:
                current_nni['runtime'] = float(line.split(":")[1].strip())
            elif "Optimal score for this search:" in line:
                optimal_score = float(line.split(":")[1].strip())

    if current_nni:
        nni_data.append(current_nni)

    for nni in nni_data:
        total_runtime += nni['runtime']
        total_trees += nni['trees']

    avg_runtime_per_tree = total_runtime / total_trees if total_trees > 0 else 0
    avg_trees_checked = statistics.mean([nni['trees'] for nni in nni_data])

    return {
        'avg_runtime_per_tree': avg_runtime_per_tree,
        'avg_trees_checked': avg_trees_checked,
        'optimal_score': optimal_score
    }

def process_log_files(directory):
    results = {}
    for filename in os.listdir(directory):
        if filename.endswith(".log"):
            file_path = os.path.join(directory, filename)
            results[filename] = parse_log_file(file_path)
    return results

# Usage
import glob
import os

log_files = glob.glob('/Users/gc3045/scmail_v1/LAML/compare_EMsolvers/*.log')

for filename in log_files:
    data = parse_log_file(filename)
    print(f"File: {filename}")
    print(f"Average runtime per tree: {data['avg_runtime_per_tree']:.6f} seconds")
    print(f"Average number of trees checked: {data['avg_trees_checked']:.2f}")
    print(f"Optimal score: {data['optimal_score']}")
    print()

