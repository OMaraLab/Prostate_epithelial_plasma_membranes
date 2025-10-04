import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import pearsonr

def read_and_flatten_file(file_path):
    data = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('@') or line.startswith('&&'):
                    continue
                data.extend(map(float, line.split()))
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None
    return data

def calculate_pearson_r(data1, data2):
    if len(data1) > 0 and len(data2) > 0:
        return pearsonr(data1, data2)[0]
    else:
        return None

def process_cell_type_files(cell_type, replicate):
    replicate_dir = f"./{cell_type}_r{replicate}"
    lipid_types = ["Sat", "Mono", "Poly"]
    lipid_data_upper = {}
    lipid_data_lower = {}
    
    # Read data for upper and lower (only for Sat, Mono, Poly)
    for lipid in lipid_types:
        file_path_upper = f"{replicate_dir}/Nor_{cell_type}_r{replicate}_{lipid}_upper.dat"
        file_path_lower = f"{replicate_dir}/Nor_{cell_type}_r{replicate}_{lipid}_lower.dat"
        
        lipid_data_upper[lipid] = read_and_flatten_file(file_path_upper)
        lipid_data_lower[lipid] = read_and_flatten_file(file_path_lower)

    if any(data is None for data in lipid_data_upper.values()) or any(data is None for data in lipid_data_lower.values()):
        print(f"Skipping replicate {replicate} for cell type {cell_type} due to missing data.")
        return None
    
    # Initialize the correlation matrix for Sat, Mono, Poly
    correlation_matrix = np.full((len(lipid_types), len(lipid_types)), np.nan)
    
    # Calculate Pearson's R for the available data
    for i, lipid1 in enumerate(lipid_types):
        for j, lipid2 in enumerate(lipid_types):
            if i <= j:
                if lipid_data_upper[lipid1] is None or lipid_data_lower[lipid2] is None:
                    correlation_matrix[i, j] = np.nan
                    correlation_matrix[j, i] = np.nan
                else:
                    R = calculate_pearson_r(lipid_data_upper[lipid1], lipid_data_lower[lipid2])
                    if R is not None:
                        correlation_matrix[i, j] = R
                        correlation_matrix[j, i] = R
                    else:
                        correlation_matrix[i, j] = np.nan
                        correlation_matrix[j, i] = np.nan
    return correlation_matrix

def plot_average_r_heatmap(cell_type, correlation_matrices, output_filename):
    # If no valid matrices exist, skip plotting
    if correlation_matrices is None or len(correlation_matrices) == 0:
        print(f"No valid correlation matrices for {cell_type}. Skipping plot.")
        return

    # Calculate the average R for each lipid pair
    average_matrix = np.nanmean(correlation_matrices, axis=0)

    labels = ["Sat", "Mono", "Poly"]

    # Define the custom colormap
    colors = ['#C23F42', '#FF9999', '#FFFFFF', '#66B2FF', '#3366FF']
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", colors, N=256)
    
    # Plot heatmap
    plt.figure(figsize=(1.77, 1.77))
    ax = sns.heatmap(
        average_matrix, annot=True, fmt=".2f", cmap=cmap,
        xticklabels=labels, yticklabels=labels,
        center=0, vmin=-1, vmax=1, cbar=True,
        linewidths=0.5, linecolor="black",
        cbar_kws={"shrink": 0.75},
        mask=np.isnan(average_matrix)
    )
    
    # Set the title
    plt.title(f"{cell_type}", fontsize=8, loc='center')
    ax.set_xticklabels(labels, fontsize=8, rotation=90)
    ax.set_yticklabels(labels, fontsize=8, rotation=0)
    
    # Label the axes for "Upper Leaflet" and "Lower Leaflet"
    ax.set_xlabel('Upper Leaflet', fontsize=8)
    ax.set_ylabel('Lower Leaflet', fontsize=8)

    # Set font size for the annotations
    for text in ax.texts:
        text.set_fontsize(8)
    
    # Configure colorbar
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([-1.0, 0.0, 1.0])
    colorbar.set_ticklabels(["-1.00", "0.00", "1.00"])
    colorbar.ax.tick_params(labelsize=8, length=0)
    ax.set_aspect("equal")
    
    # Save the figure
    plt.savefig(output_filename, dpi=600, bbox_inches="tight")
    plt.close()

def process_all_files():
    cell_types = ["RWPE-1", "PS_PC-3", "PS_DU145", "PS_LNCaP", "BPH-1", "PC-3", "DU145", "LNCaP"]
    replicates = [1, 2, 3]
    
    for cell_type in cell_types:
        correlation_matrices = []
        for replicate in replicates:
            correlation_matrix = process_cell_type_files(cell_type, replicate)
            if correlation_matrix is not None:
                correlation_matrices.append(correlation_matrix)
        
        if correlation_matrices:
            correlation_matrices = np.array(correlation_matrices)
            
            # Plot and save the average R heatmap
            output_png = f"heatmap_{cell_type}_upperLower.png"
            plot_average_r_heatmap(cell_type, correlation_matrices, output_png)
            print(f"Saved heatmap for {cell_type} to {output_png}")
        else:
            print(f"No valid data for cell type {cell_type}, skipping heatmap generation.")

process_all_files()

