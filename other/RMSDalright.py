from vina import Vina 
from datetime import datetime
import numpy as np
import pandas as pd
import subprocess
import os

# Function to calculate RMSD using OpenBabel
def calculate_rmsd(reference, ligand, n_poses=10):
    ref_mol2 = reference.replace(".pdbqt", ".mol2")
    lig_mol2 = ligand.replace(".pdbqt", ".mol2")

    # Convert reference and ligand to MOL2 format with hydrogens
    subprocess.run(f"obabel {reference} -O {ref_mol2} --gen3D --addh", shell=True, check=False)
    subprocess.run(f"obabel {ligand} -O {lig_mol2} --gen3D --addh", shell=True, check=False)

    rmsd_values = []
    try:
        # Compute RMSD using OpenBabel for all poses in the single file
        result = subprocess.run(f"obrms -f {ref_mol2} {lig_mol2}", shell=True, capture_output=True, text=True, check=True)

        # Extract RMSD values from output
        output_lines = result.stdout.strip().split("\n")
        for line in output_lines[:n_poses]:  # Limit to expected number of poses
            parts = line.split()
            if parts and parts[-1].replace('.', '', 1).isdigit():  # Ensure last part is a number
                rmsd_values.append(float(parts[-1]))

        #Ignore unreasonable RMSD values
        rmsd_values = [r for r in rmsd_values if r < 50]

    except Exception as e:
        print(f" Error computing RMSD for {ligand}: {e}")

    # Clean up temporary files
    os.remove(ref_mol2)
    os.remove(lig_mol2)

    return rmsd_values

# Define docking function
def dock_single_point(center):
    v = Vina(sf_name='vina')
    v.set_receptor('8IBT_receptor.pdbqt')
    v.set_ligand_from_file('LNT_ligand.pdbqt')

    box_size = [18, 10, 13]
    affinities = []
    rmsd_values = []

    for repeat in range(1):  # Run 3 docking repeats
        v.compute_vina_maps(center=center, box_size=box_size)
        v.dock(exhaustiveness=32, n_poses=10)

        output_file = f'ligand_docked_{center[0]}_{center[1]}_{center[2]}_rep{repeat + 1}.pdbqt'
        v.write_poses(output_file, n_poses=10, overwrite=True)

        energy = v.energies()
        if energy.size > 0:
            affinities.append(energy[0, 0])

    # Compute RMSD for all poses in the single docking file
    rmsd_values = calculate_rmsd('LNT_ligand.pdbqt', output_file)

    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    if affinities:
        mean_aff = np.mean(affinities)
        std_aff = np.std(affinities)
        aff_str = ', '.join(f"{a:.2f}" for a in affinities)

        # Store RMSD values as a single string
        rmsd_str = ', '.join(f"{r:.2f}" for r in rmsd_values) if rmsd_values else 'N/A'

        result = [timestamp, str(center), mean_aff, std_aff, aff_str, rmsd_str]
    else:
        result = [timestamp, str(center), 'No Pose', 'N/A', 'N/A', 'N/A']

    return result

# Run docking for the single coordinate (75, -30, -60)
grid_center = (75, -30, -60)

print(f"Starting docking for grid center: {grid_center}")
start_time = datetime.now()
result = dock_single_point(grid_center)
end_time = datetime.now()
print(f"Execution time: {end_time - start_time}")

# Save results to file
with open('affinity_results_single.txt', 'w') as result_file:
    result_file.write('Timestamp\tGrid Center\tMean Affinity\tStd Dev\tAll Affinities\tRMSDs per Pose\n')
    result_file.write('---------\t-----------\t--------------\t--------\t----------------\t--------------\n')
    result_file.write('\t'.join(map(str, result)) + '\n')

df = pd.DataFrame([result], columns=["Timestamp", "Grid Center", "Mean Affinity", "Std Dev", "Affinities", "RMSDs"])
df.to_csv("affinity_results_single.csv", index=False)

print("\n?? Docking for single coordinate completed!")
print("Results saved in 'affinity_results_single.txt' and 'affinity_results_single.csv'")
