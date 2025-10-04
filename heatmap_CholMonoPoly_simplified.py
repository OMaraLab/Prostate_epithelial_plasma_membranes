import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import matplotlib.colors as mcolors

def read_correlation_matrix(filename):
    """ Reads a 3x3 correlation matrix from a .dat file, skipping lines with '@' or '&&'. """
    with open(filename, 'r') as file:
        lines = [line.strip() for line in file if not line.startswith(("@", "&&"))]
    
    # Convert lines to a NumPy array
    data = np.array([list(map(float, line.split())) for line in lines])
    return data

def plot_correlation_heatmap(correlation_matrix, output_filename):
    """ Generates and saves a heatmap for the given 3x3 correlation matrix. """
    labels = ["Mono", "Poly", "CHOL"]

    # Define the custom colormap with 5 colours (transition at -0.75, -0.5, 0.0, 0.5, 0.75)
    colors = ['#0000FF', '#99CCFF', '#FFFFFF', '#FF9999', '#FF0000']  # Primary blue, pastel blue, white, pastel red, bright red
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", colors, N=256)

    # Create heatmap
    plt.figure(figsize=(1.97, 2.36))  # 5 cm × 6 cm in inches
    ax = sns.heatmap(
        correlation_matrix, annot=True, fmt=".2f",
        xticklabels=labels, yticklabels=labels,
        cmap=cmap, center=0, vmin=-1, vmax=1, cbar=True,
        linewidths=0.5, linecolor="black",
        cbar_kws={"shrink": 0.75}  # Adjust color bar size
    )

    # Set font size and rotate labels
    ax.set_xticklabels(labels, fontsize=8, rotation=90)  # Column labels rotated vertically
    ax.set_yticklabels(labels, fontsize=8, rotation=0)   # Row labels rotated horizontally
    for text in ax.texts:  
        text.set_fontsize(8)  # Set annotation font size

    # Color bar customization
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([-1.0, 0.0, 1.0])  # Only display ticks at -1, 0, and 1
    colorbar.set_ticklabels(["-1.00", "0.00", "1.00"])  # Ticks with two decimal places
    colorbar.ax.tick_params(labelsize=8, length=0)  # Remove tick marks by setting length=0

    # Ensure square aspect ratio
    ax.set_aspect("equal")

    # Save the plot
    plt.savefig(output_filename, dpi=600, bbox_inches="tight")
    plt.close()

# Process all correlation matrix files
files = glob.glob("correlation_*_r1.dat")  # Update to match the new filenames
for file in files:
    matrix = read_correlation_matrix(file)
    output_png = file.replace(".dat", ".png")  # Save as PNG with the same filename
    plot_correlation_heatmap(matrix, output_png)
    print(f"Saved: {output_png}")

