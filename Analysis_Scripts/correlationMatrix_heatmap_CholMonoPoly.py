import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import matplotlib.colors as mcolors
from scipy.stats import pearsonr

# Function to read and flatten the data from a file
def read_and_flatten_file(file_path):
    data = []
    try:
        print(f"Processing file: {file_path}")
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('@') or line.startswith('&&'):
                    continue
                data.extend(map(float, line.split()))
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None
    return data

# Function to calculate the Pearson correlation coefficient
def calculate_pearson_r(data1, data2):
    if len(data1) != len(data2):
        print("Data arrays must have the same length!")
        return None
    return pearsonr(data1, data2)[0]

# Function to process each cell type's files and compute correlation matrix
def process_cell_type_files(cell_type, replicate):
    replicate_dir = f"./{cell_type}_r{replicate}"
    
    lipid_types = ["CHOL", "Mono", "Poly"]
    lipid_data = {}

    for lipid in lipid_types:
        file_path = f"{replicate_dir}/Nor_{cell_type}_r{replicate}_{lipid}.dat"
        lipid_data[lipid] = read_and_flatten_file(file_path)
    
    if any(data is None for data in lipid_data.values()):
        print(f"Skipping {cell_type}_r{replicate} due to missing data.")
        return None

    correlation_matrix = np.ones((len(lipid_types), len(lipid_types)))  # Diagonal elements are 1

    for i, lipid1 in enumerate(lipid_types):
        for j, lipid2 in enumerate(lipid_types):
            if i < j:  # Only compute upper triangle
                R = calculate_pearson_r(lipid_data[lipid1], lipid_data[lipid2])
                if R is not None:
                    correlation_matrix[i, j] = R
                    correlation_matrix[j, i] = R  # Symmetric matrix

    return correlation_matrix

# Function to plot the heatmap for a cell type based on average correlation
def plot_correlation_heatmap(cell_type, correlation_matrices, output_filename):
    mean_matrix = np.mean(correlation_matrices, axis=0)
    std_matrix = np.std(correlation_matrices, axis=0)

    labels = ["Chol", "Mono", "Poly"]

    colors = ['#C23F42', '#FF9999', '#FFFFFF', '#66B2FF', '#3366FF']
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", colors, N=256)

    plt.figure(figsize=(1.5, 1.5))  
    ax = sns.heatmap(
        mean_matrix, annot=True, fmt=".2f", cmap=cmap,
        xticklabels=labels, yticklabels=labels,
        center=0, vmin=-1, vmax=1, cbar=True,
        linewidths=0.5, linecolor="black",
        cbar_kws={"shrink": 0.75}
    )

    plt.title(cell_type, fontsize=7, loc='center')

    ax.set_xticklabels(labels, fontsize=7, rotation=90)
    ax.set_yticklabels(labels, fontsize=7, rotation=0)
    for text in ax.texts:
        text.set_fontsize(7)

    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([-1.0, 0.0, 1.0])
    colorbar.set_ticklabels(["-1.00", "0.00", "1.00"])
    colorbar.ax.tick_params(labelsize=7, length=0)

    ax.set_aspect("equal")

    plt.savefig(output_filename, dpi=600, bbox_inches="tight")
    plt.close()

# Function to process all replicates and generate the heatmap for each cell type
def process_all_files():
    cell_types = ["RWPE-1", "PS_PC-3", "PS_DU145", "PS_LNCaP", "BPH-1", "PC-3", "DU145", "LNCaP"]
    replicates = [1, 2, 3]
    
    for cell_type in cell_types:
        correlation_matrices = []
        
        for replicate in replicates:
            correlation_matrix = process_cell_type_files(cell_type, replicate)
            if correlation_matrix is not None:
                correlation_matrices.append(correlation_matrix)
        
        if len(correlation_matrices) > 0:
            correlation_matrices = np.array(correlation_matrices)
            output_png = f"heatmap_{cell_type}.png"
            plot_correlation_heatmap(cell_type, correlation_matrices, output_png)
            print(f"Saved heatmap for {cell_type} to {output_png}")

# Run the processing
process_all_files()

