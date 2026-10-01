import os
import re
import numpy as np
import matplotlib.pyplot as plt

# --------- SETTINGS ---------

stiffnesses = [0.5, 1, 2, 4, 8, 10, 12, 14, 16, 32]
N_REP = 3  # number of replicates

# --------- FILE HANDLING ---------

def get_last_timestep(directory):
    max_step = -1

    for fname in os.listdir(directory):
        match = re.match(r"ctags(\d+)_r\d+\.out$", fname)
        if match:
            step = int(match.group(1))
            if step > max_step:
                max_step = step

    if max_step < 0:
        raise ValueError(f"No MPI ctags files in {directory}")

    return max_step

def load_ctags_mpi(directory, step):
    import glob

    pattern = os.path.join(directory, f"ctags{step}_r*.out")
    files = glob.glob(pattern)

    if not files:
        raise ValueError(f"No files for step {step} in {directory}")

    # sort by rank
    def get_rank(f):
        m = re.search(r"_r(\d+)", f)
        return int(m.group(1)) if m else 0

    files = sorted(files, key=get_rank)

    chunks = []
    for f in files:
        data = np.loadtxt(f)

        if data.ndim == 1:
            data = data.reshape(1, -1)

        chunks.append(data)

    # stack in y direction
    full = np.hstack(chunks)

    # convert to (Y, X)
    full = full.T

    return full

# --------- DATA PROCESSING ---------

def load_points(filepath):
    grid = np.loadtxt(filepath)
    return np.argwhere(grid > 0)


def compute_area(points):
    return len(points)


def compute_length(points):
    max_dist = 0.0
    n = len(points)

    for i in range(n):
        for j in range(i + 1, n):
            dx = points[i, 0] - points[j, 0]
            dy = points[i, 1] - points[j, 1]
            dist = np.sqrt(dx * dx + dy * dy)

            if dist > max_dist:
                max_dist = dist

    return max_dist


def measure_directory(directory):
    step = get_last_timestep(directory)
    grid = load_ctags_mpi(directory, step)
    
    # sanity check the grid space
    #print(f"  grid shape: {grid.shape}")

    points = np.argwhere(grid > 0)

    area = compute_area(points)
    length = compute_length(points)

    return area, length


# --------- MAIN ---------

mean_areas = []
mean_lengths = []

std_areas = []
std_lengths = []

print("\nProcessing simulations...\n")

for stiffness in stiffnesses:
    areas = []
    lengths = []

    print(f"Stiffness {stiffness} kPa")

    for rep in range(1, N_REP + 1):
        d = f"../output/run_{stiffness}_kPa_rep{rep}"

        if not os.path.isdir(d):
            print(f"  WARNING: missing {d}")
            continue

        area, length = measure_directory(d)

        areas.append(area)
        lengths.append(length)

        print(f"  rep{rep}: area={area}, length={length}")

    if len(areas) == 0:
        mean_areas.append(np.nan)
        mean_lengths.append(np.nan)
        std_areas.append(np.nan)
        std_lengths.append(np.nan)
        continue

    mean_areas.append(np.mean(areas))
    mean_lengths.append(np.mean(lengths))

    std_areas.append(np.std(areas))
    std_lengths.append(np.std(lengths))


# --------- PLOTS ---------

# Area
plt.figure()
positions = np.arange(len(stiffnesses))
plt.errorbar(positions, mean_areas, yerr=std_areas,
             marker='o', color='black', linestyle='none', capsize=5)
plt.xticks(positions, labels=[str(s) for s in stiffnesses])
plt.xlabel("Stiffness (kPa)")
plt.ylabel("Cell Area (pixels)")
plt.title("Cell Area vs Stiffness")
plt.grid(False)
plt.savefig("cell_area_vs_stiffness.png")

# Length
plt.figure()
positions = np.arange(len(stiffnesses))
plt.errorbar(positions, mean_lengths, yerr=std_lengths,
             marker='o', color='black', linestyle='none', capsize=5)
plt.xticks(positions, labels=[str(s) for s in stiffnesses])
plt.xlabel("Stiffness (kPa)")
plt.ylabel("Cell Length (pixels)")
plt.title("Cell Length vs Stiffness")
plt.grid(False)
plt.savefig("cell_length_vs_stiffness.png")