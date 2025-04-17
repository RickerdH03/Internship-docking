from vina import Vina 
from datetime import datetime
import numpy as np
import pandas as pd
import multiprocessing
from tqdm import tqdm

# Generate exact axis values (5 steps, includes start and end)
def linspace_axis(start, end, num_points=5):
    return [round(v, 2) for v in np.linspace(start, end, num_points)]

# Define docking function for multiprocessing
def dock_grid_point(center):
    v = Vina(sf_name='vina')
    v.set_receptor('8IBT_receptor.pdbqt')
    v.set_ligand_from_file('LNT_ligand.pdbqt')

    box_size = [18, 10, 13]
    affinities = []

    for repeat in range(3):
        v.compute_vina_maps(center=center, box_size=box_size)
        v.dock(exhaustiveness=8, n_poses=1)

        output_file = f'ligand_docked_{center[0]}_{center[1]}_{center[2]}_rep{repeat + 1}.pdbqt'
        v.write_poses(output_file, n_poses=1, overwrite=True)

        energy = v.energies()
        if energy.size > 0:
            affinities.append(energy[0, 0])

    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    if affinities:
        mean_aff = np.mean(affinities)
        std_aff = np.std(affinities)
        aff_str = ', '.join(f"{a:.2f}" for a in affinities)
        result = [timestamp, str(center), mean_aff, std_aff, aff_str]
    else:
        result = [timestamp, str(center), 'No Pose', 'N/A', 'N/A']

    return result

# Define axis ranges
x_vals = linspace_axis(60.36, 79.24, 5)
y_vals = linspace_axis(-37.89, -14.18, 5)
z_vals = linspace_axis(-71.77, -52.21, 5)

# Generate all grid centers
grid_centers = [(x, y, z) for x in x_vals for y in y_vals for z in z_vals]

# Sequential docking
#for center in tqdm(grid_centers):
#    dock_grid_point(center)

# Parallel docking
num_workers = multiprocessing.cpu_count()
print(f"Number of workers (CPUs): {num_workers}")

start_time = datetime.now()
with multiprocessing.Pool(processes=num_workers) as pool:
    results_parallel = list(tqdm(pool.imap(dock_grid_point, grid_centers), total=len(grid_centers)))
end_time = datetime.now()
print(f"Parallel execution time: {end_time - start_time}")

 # Save results to files
    with open('affinity_results.txt', 'w') as result_file:
        result_file.write('Timestamp\tGrid Center\tMean Affinity\tStd Dev\tAll Affinities\n')
        result_file.write('---------\t-----------\t--------------\t--------\t----------------\n')
        for row in results:
            result_file.write('\t'.join(map(str, row)) + '\n')

    df = pd.DataFrame(results, columns=["Timestamp", "Grid Center", "Mean Affinity", "Std Dev", "Affinities"])
    df.to_csv("affinity_results.csv", index=False)

    print("\nðŸŽ¯ All docking jobs completed in parallel!")
    print("Results saved in 'affinity_results.txt' and 'affinity_results.csv'")