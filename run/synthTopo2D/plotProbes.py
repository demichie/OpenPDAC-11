import matplotlib.pyplot as plt
import numpy as np

# Function to read the text file and extract data
def read_text_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        data = [line.strip().split() for line in lines if not line.startswith("#")]
    return data

# Function to plot the data
def plot_data(file_path):
    data = read_text_file(file_path)

    times = [float(row[0]) for row in data]
    values = [list(map(float, row[1:])) for row in data]

    num_probes = len(values[0])
    for i in range(num_probes):
        plt.plot(times, [row[i] for row in values], label=f"Probe {i + 1}")

    plt.xlabel('Time (s)')
    plt.ylabel(file_path)
    plt.legend()
    plt.show()

# Replace 'file_path' with your actual file path
file_path = 'postProcessing/probes/0/p'
plot_data(file_path)

