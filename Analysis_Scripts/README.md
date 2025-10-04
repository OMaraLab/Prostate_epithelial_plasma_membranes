# Code Usage Guide

## Example Workflow

This workflow demonstrates a complete analysis pipeline from trajectory analysis to density correlation:

1. **Analyse membrane properties** (Section 1-2): Calculate thickness, APL, tilt angles, flip-flop events, and diffusion coefficients from trajectory files
2. **Generate density maps** (Section 3): Use g_mydensity to create lipid number density heatmaps
3. **Visualise density data** (Section 4): Create heatmap images with dispgrid
4. **Calculate correlations** (Sections 5-7): Compute spatial co-localisation between lipid types
5. **Generate final heatmaps** (Section 8): Create correlation heatmaps

## Membrane Property Analysis

### 1. compute_membrane_properties.ipynb, for calculating structural membrane properties using Lipyds

Computes four membrane properties using the Lipyds package (https://github.com/lilyminium/2023-10-10_lipyd-example) with all data saved to `.pkl` files:

**1. Membrane thickness** - The `MembraneThickness` class calculates thickness as mean distance between upper and lower leaflet phosphate groups over a 2D grid. The code computes grid-based averages and reports overall mean thickness.

**2. Area per lipid (APL)** - The `AreaPerLipid` class uses `scipy.spatial.Voronoi` to compute 2D Voronoi tessellation:

$$
\text{APL}_i = A_{\text{Voronoi},i}
$$

where $A_{\text{Voronoi},i}$ is the polygon area surrounding lipid $i$ after projecting neighbours onto the local membrane plane. Code filters values >200 Å².

**3. Cholesterol tilt angle** - The `LipidTilt` class calculates the cholesterol tilt angles.

**4. Cholesterol flip-flop events** - The `LipidFlipFlop` class assigns cholesterol to leaflets based on minimum distance to headgroups. If minimum distance < `leaflet_width`, cholesterol is in that leaflet; otherwise interstitial. Events counted for direct upper↔lower transitions.

**Example usage:**
```python
import MDAnalysis as mda
from lipyds.leafletfinder.leafletfinder import LeafletFinder
from lipyds.analysis.thickness import MembraneThickness
from lipyds import AreaPerLipid
from lipyds.analysis.tilt import LipidTilt
from lipyds.analysis.flipflop import LipidFlipFlop

u = mda.Universe("topol.tpr", "traj.xtc")
membrane_residues = u.select_atoms("not resname ION W WF")

finder = LeafletFinder(universe=membrane_residues, select="name PO4 GL1 GL2 AM1 AM2",
                       cutoff=30, method="spectralclustering", n_leaflets=2)

# Thickness
thickness = MembraneThickness(universe=membrane_residues, select="name PO4",
                              leafletfinder=finder, bin_size=4)
thickness.run()

# APL
apl = AreaPerLipid(universe=membrane_residues, select="name PO4 ROH ES",
                   leafletfinder=finder, cutoff=30)
apl.run()

# Tilt
chol = u.select_atoms("resname CHOL")
tilt = LipidTilt(universe=chol, leafletfinder=finder,
                 select="name R1", select_end="name R5", cutoff=10)
tilt.run()

# Flip-flop
flipflop = LipidFlipFlop(universe=chol, leafletfinder=finder,
                         select="name R2", leaflet_width=11, cutoff=5)
flipflop.run()
```

**Directory structure:**
```
project_root/
├── compute_membrane_properties.ipynb
├── topol.tpr
└── traj.xtc
```

**Output:** `THICKNESS.pkl`, `RESNAME_APL_df.pkl`, `Tilt_result_df.pkl`, `flipflop_results.pkl`, `resname_mean_std.xlsx`, `LF_resname_mean_std.xlsx`.

### 2. compute_membrane_lateral_diffusion.ipynb, for calculating lipid diffusion coefficients using Lipyphilic

Calculates lateral diffusion coefficients via MSD analysis using the Lipyphilic package (https://pubs.acs.org/doi/10.1021/acs.jctc.1c00447). Applies periodic boundary unwrapping, computes MSD vs. lag time, and extracts diffusion coefficients by linear regression.

**MSD calculation:**

The `MSD` class computes:

$$
\text{MSD}_i(\tau) = \frac{1}{N_{\text{origins}}} \sum_{t} [x_i(t+\tau) - x_i(t)]^2 + [y_i(t+\tau) - y_i(t)]^2
$$

where $x_i(t)$ and $y_i(t)$ are membrane plane coordinates of lipid $i$ at time $t$, $\tau$ is lag time, summed over all valid time origins.

**Diffusion coefficient extraction:**

The `diffusion_coefficient` method uses `scipy.stats.linregress` to fit the linear regime (100–1000 ns):

$$
\text{MSD}(\tau) = 4D\tau + c
$$

$$
D = \frac{1}{4} \times \text{slope}
$$

where the factor of 4 arises from 2D diffusion, and slope is extracted from linear regression of MSD vs. $\tau$.

**Example usage:**
```python
import MDAnalysis as mda
from lipyphilic.transformations import nojump
from lipyphilic.lib.lateral_diffusion import MSD

u = mda.Universe("topol.tpr", "traj.xtc")
ag = u.select_atoms("name GL1 GL2 AM1 AM2 ROH ES")

u.trajectory.add_transformations(nojump(ag=ag, nojump_x=True, nojump_y=True, nojump_z=False))

msd = MSD(universe=u, lipid_sel="name GL1 GL2 AM1 AM2 ES ROH",
          com_removal_sel="name GL1 GL2 AM1 AM2 ES ROH")
msd.run()

# Whole membrane
d, sem = msd.diffusion_coefficient(start_fit=100, stop_fit=1000)

# Specific selections
d_upper, sem_upper = msd.diffusion_coefficient(start_fit=100, stop_fit=1000,
                                               lipid_sel="name GL1 GL2 AM1 AM2 and resid 1 to 1900")

d_pc, sem_pc = msd.diffusion_coefficient(start_fit=100, stop_fit=1000, lipid_sel="resname *PC")
```

**Directory structure:**
```
project_root/
├── compute_membrane_lateral_diffusion.ipynb
├── topol.tpr
└── traj.xtc
```

**Output:** `MSD.pkl`, `Membrane_UL_MSD.xlsx`, `Membrane_LL_MSD.xlsx`.

Analyses whole membrane, individual leaflets (by residue ID), lipid species (wildcards like `*PC`), and saturation classes. For multiple replicates, average MSD arrays before calculating diffusion coefficients.

## Density Map Generation and Visualisation

### 3. g_mydensity (GROMACS 4.5.7), for generating lipid number density data from molecular dynamics trajectory files

Use g_mydensity (https://doi.org/10.1016/j.chemphyslip.2013.02.001) to generate lipid number density data from simulation trajectories. First, create an index file grouping lipid types to analyse. The output `.dat` file contains raw number density values that must be normalised for comparison between systems:

$$
\text{Normalised density}(x,y) = \frac{\text{Number density}(x,y)}{N_{\text{lipid}}}
$$

where $\text{Number density}(x,y)$ is the count of lipids in heatmap bin $(x,y)$, and $N_{\text{lipid}}$ is the total number of that lipid species in the membrane or leaflet.

The $N_{\text{lipid}}$ information is in the **factors** folder

**Example usage:**
```bash
gmx make_ndx -f conf.gro -o lipids.ndx
g_mydensity -f traj.xtc -s topol.tpr -n lipids.ndx -o CHOL.xvg -og CHOL.dat -dens number
# Manually normalise and save as Nor_Membrane1_r1_CHOL.dat
```

**Directory structure:**
```
project_root/
├── Membrane1_r1/dat_files/
│   ├── Membrane1_r1_CHOL.dat      # Raw output
│   └── ...
└── factors.txt                 # N_lipid for each type
```

### 4. dispgrid_blue_white_yellow_red.py, for heatmap visualisation from normalised density data

Converts normalised density grid data into heatmap images using a custom blue-white-yellow-red colourmap. Accepts data from g_mydensity or similar tools (g_thickness). Based on dispgrid.py from https://doi.org/10.1016/j.chemphyslip.2013.02.001.

**Example usage:**
```bash
python dispgrid_blue_white_yellow_red.py Nor_BPH-1_r1_CHOL.dat CHOL_heatmap.png
python dispgrid_blue_white_yellow_red.py Nor_BPH-1_r1_CHOL.dat CHOL_heatmap.png \
    --zlim 0 1.0 --nlevels 3
```

Key parameters: `--sampling` (sampling file), `--min_fraction` (default 0.5), `--zlim` (colour scale limits), `--xlim/--ylim` (axis limits), `--format` (output format), `--nlevels` (colour levels), `--discrete` (discrete vs interpolated), `--center` (centre origin).

## Spatial Correlation Analysis

### 5. correlationMatrix_heatmap_all.py, for comprehensive spatial co-localisation analysis across all lipid types

Calculates Pearson correlation coefficients between lipid density distributions to identify spatial co-localisation patterns. Analyses nine lipid categories: CHOL, Sat, Mono, Poly, PC, PE, SM, PS, PI.

The code uses `scipy.stats.pearsonr` which implements:

$$
R_{ij} = \frac{\sum_{k=1}^{n} (d_{i,k} - \bar{d_i})(d_{j,k} - \bar{d_j})}{\sqrt{\sum_{k=1}^{n}(d_{i,k} - \bar{d_i})^2} \sqrt{\sum_{k=1}^{n}(d_{j,k} - \bar{d_j})^2}}
$$

where $d_{i,k}$ is normalised density of lipid type $i$ in bin $k$, $\bar{d_i} = \frac{1}{n}\sum_{k=1}^{n}d_{i,k}$ is the mean density, and $n$ is total bins. Range: -1 (perfect negative) to +1 (perfect positive correlation).

**Example usage:**
```bash
python correlationMatrix_heatmap_all.py  # Auto-processes all cell types/replicates
```

**Directory structure:**
```
project_root/
├── correlationMatrix_heatmap_all.py
├── BPH-1_r1/
│   ├── Nor_BPH-1_r1_CHOL.dat
│   ├── Nor_BPH-1_r1_Sat.dat
│   └── ... (all 9 lipid types)
├── BPH-1_r2/ [same]
└── BPH-1_r3/ [same]
```

**Output:** Individual replicate matrices (`correlation_matrix_BPH-1_r{1,2,3}.txt`), averaged matrix (`correlation_matrix_BPH-1_all.txt`), heatmap (`heatmap_BPH-1_all.png`).

### 6. correlationMatrix_heatmap_CholMonoPoly.py, for generation of lipid density correlation matrices where the Pearson's R values are for tail classifications

Subset analysis of CHOL, Mono, and Poly only using `scipy.stats.pearsonr`. Identical methodology to script 5 but produces 3×3 correlation heatmaps averaged across replicates.

**Directory structure:** Same as script 5 but requires only `Nor_{cell}_r{rep}_{CHOL|Mono|Poly}.dat` files.

**Output:** `heatmap_BPH-1.png` (3×3 matrix).

### 7. correlationMatrix_heatmap_upperLower.py, for transbilayer lipid co-localisation between leaflets (interdigitation)

Analyses lipid registration across bilayer by correlating densities between upper and lower leaflets using `scipy.stats.pearsonr`. Examines transbilayer coupling rather than within-leaflet co-localisation.

$$
R_{ij}^{\text{trans}} = \frac{\sum_{k=1}^{n} (d_{i,k}^{\text{upper}} - \bar{d_i^{\text{upper}}})(d_{j,k}^{\text{lower}} - \bar{d_j^{\text{lower}}})}{\sqrt{\sum_{k=1}^{n}(d_{i,k}^{\text{upper}} - \bar{d_i^{\text{upper}}})^2} \sqrt{\sum_{k=1}^{n}(d_{j,k}^{\text{lower}} - \bar{d_j^{\text{lower}}})^2}}
$$

where superscripts indicate leaflet and $k$ indexes corresponding bins on opposite leaflets.

**Directory structure:** Requires separate `Nor_{cell}_r{rep}_{Sat|Mono|Poly}_{upper|lower}.dat` files.

**Output:** `heatmap_BPH-1_upperLower.png` (cross-leaflet correlations).

### 8. heatmap_CholMonoPoly_simplified.py, for rapid visualisation of pre-calculated matrices

Generates heatmaps from existing 3×3 correlation matrix files without recalculation. Reads `correlation_*_r1.dat` files (skipping lines with `@` or `&&`).

**Input format:**
```
@ metadata (skipped)
1.00  -0.45  0.23
-0.45  1.00  0.67
0.23  0.67  1.00
```

**Output:** `correlation_*_r1.png`.
