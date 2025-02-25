import re
import argparse

def parse_log_file(filename):
    oldEM_data = []
    fastEM_data = []

    pattern = r"Current (\w+)_solver score: ([-\d.]+), ([-\d.]+), ([-\d.]+), (\d+\.\d+) seconds"

    with open(filename, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                solver, llh, phi, nu, time = match.groups()
                data = (float(llh), float(phi), float(nu), float(time))
                if solver == 'oldEM':
                    oldEM_data.append(data)
                elif solver == 'fastEM':
                    fastEM_data.append(data)

    return oldEM_data, fastEM_data

def calculate_statistics(oldEM_data, fastEM_data):
    differences = []
    for old, fast in zip(oldEM_data, fastEM_data):
        diff = [old[i] - fast[i] for i in range(3)]  # Only llh, phi, nu
        differences.append(diff)

    avg_differences = [sum(col) / len(col) for col in zip(*differences)]
    
    oldEM_avg_time = sum(data[3] for data in oldEM_data) / len(oldEM_data)
    fastEM_avg_time = sum(data[3] for data in fastEM_data) / len(fastEM_data)

    return avg_differences, oldEM_avg_time, fastEM_avg_time

def main():
    parser = argparse.ArgumentParser(description="Analyze EM solver log file")
    parser.add_argument("logfile", help="Path to the log file")
    args = parser.parse_args()

    oldEM_data, fastEM_data = parse_log_file(args.logfile)
    avg_differences, oldEM_avg_time, fastEM_avg_time = calculate_statistics(oldEM_data, fastEM_data)

    print(f"Average differences (oldEM - fastEM):")
    print(f"LLH: {avg_differences[0]:.6f}")
    print(f"Phi: {avg_differences[1]:.6f}")
    print(f"Nu: {avg_differences[2]:.6f}")
    print(f"Average time for oldEM_solver: {oldEM_avg_time:.6f} seconds")
    print(f"Average time for fastEM_solver: {fastEM_avg_time:.6f} seconds")

if __name__ == "__main__":
    main()

