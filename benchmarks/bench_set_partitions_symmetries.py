import time
import matplotlib.pyplot as plt
from SBB.PyMoments.Combinatorics import set_partitions, set_partitions_symmetries

def benchmark_set_partitions():
    sizes = list(range(2, 15, 2))  # Set sizes: 2, 4, 6, 8, 10, 12
    set_templates = [['X'] * (n // 2) + ['Y'] * (n // 2) for n in sizes]

    times_set_partitions = []
    times_set_partitions_symmetries = []
    speedup_ratios = []

    for test_set in set_templates:
        tuple_input = tuple(test_set)

        # Time set_partitions
        start_time = time.time()
        list(set_partitions(tuple_input))
        end_time = time.time()
        time_partitions = end_time - start_time
        times_set_partitions.append(time_partitions)

        # Time set_partitions_symmetries
        start_time = time.time()
        set_partitions_symmetries(test_set)
        end_time = time.time()
        time_symmetries = end_time - start_time
        times_set_partitions_symmetries.append(time_symmetries)

        # Compute speedup
        speedup = time_partitions / time_symmetries if time_symmetries > 0 else float('inf')
        speedup_ratios.append(speedup)

        print(f"Set size {len(test_set)} | set_partitions: {time_partitions:.4f}s | "
              f"set_partitions_symmetries: {time_symmetries:.4f}s | Speedup: {speedup:.2f}x")

    # Plot execution times
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(sizes, times_set_partitions, label="set_partitions", marker='o')
    plt.plot(sizes, times_set_partitions_symmetries, label="set_partitions_symmetries", marker='s')
    plt.xlabel("Set Size")
    plt.ylabel("Execution Time (seconds)")
    plt.title("Execution Time: set_partitions vs. set_partitions_symmetries")
    plt.legend()
    plt.grid(True)

    # Plot speedup
    plt.subplot(1, 2, 2)
    plt.plot(sizes, speedup_ratios, label="Speedup (set_partitions / set_partitions_symmetries)", marker='^', color='red')
    plt.xlabel("Set Size")
    plt.ylabel("Speedup Factor")
    plt.title("Speedup of set_partitions_symmetries")
    plt.axhline(y=1, color='gray', linestyle='--', linewidth=0.8)  # Baseline at y=1 (no speedup)
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

def sum_multiplicities(partitions):
    """Sum the first element (count) of each partition tuple only when it's greater than 1."""
    return sum(count for count, *partition in partitions if count > 1)

def benchmark_multiplicities():
    sizes = list(range(2, 17, 2))  # Set sizes: 2, 4, 6, 8, 10, 12
    set_templates = [['X'] * (n // 2) + ['Y'] * (n // 2) for n in sizes]

    multiplicities_sums = []
    execution_times = []

    for test_set in set_templates:
        start_time = time.time()
        partitions_symmetries_result = set_partitions_symmetries(test_set)
        end_time = time.time()

        execution_time = end_time - start_time
        multiplicity_sum = sum_multiplicities(partitions_symmetries_result)

        execution_times.append(execution_time)
        multiplicities_sums.append(multiplicity_sum)

        print(f"Set size {len(test_set)} | Multiplicity Sum: {multiplicity_sum} | Time: {execution_time:.4f}s")

    # Plot multiplicity sums
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(sizes, multiplicities_sums, label="Sum of Multiplicities (>1)", marker='o', color='blue')
    plt.xlabel("Set Size")
    plt.ylabel("Sum of Multiplicities")
    plt.title("Sum of Multiplicities in set_partitions_symmetries")
    plt.legend()
    plt.grid(True)

    # Plot execution times
    plt.subplot(1, 2, 2)
    plt.plot(sizes, execution_times, label="Execution Time", marker='s', color='red')
    plt.xlabel("Set Size")
    plt.ylabel("Execution Time (seconds)")
    plt.title("Execution Time of set_partitions_symmetries")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    benchmark_set_partitions()
    benchmark_multiplicities()
