from vina import Vina
from datetime import datetime
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.PDBIO import PDBIO
import tempfile
import os
from scipy.spatial.distance import cdist

# Load the true reference ligand from an external clean PDB file
REFERENCE_LIGAND_PATH = '20250213-LNT from PDB.pdb'

def get_structure_from_atoms(atom_lines):
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w") as temp:
        temp.writelines(atom_lines)
        temp_path = temp.name

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("lig", temp_path)
    os.remove(temp_path)
    return structure

def load_structure_from_file(filepath):
    parser = PDBParser(QUIET=True)
    return parser.get_structure("ref", filepath)

def extract_heavy_atoms(structure):
    return [atom for atom in structure.get_atoms() if atom.element != 'H']

def match_atoms_by_proximity(ref_atoms, mobile_atoms):
    ref_coords = np.array([a.coord for a in ref_atoms])
    mob_coords = np.array([a.coord for a in mobile_atoms])

    dist_matrix = cdist(ref_coords, mob_coords)
    indices = dist_matrix.argmin(axis=1)

    matched_mobile = [mobile_atoms[i] for i in indices]
    return matched_mobile

def calculate_rmsd_with_superimposer(ref_structure, docked_structures):
    ref_atoms = extract_heavy_atoms(ref_structure)
    rmsds = []

    for struct in docked_structures:
        mobile_atoms = extract_heavy_atoms(struct)

        if len(ref_atoms) == 0 or len(mobile_atoms) == 0:
            print("Empty atom list. Skipping pose.")
            continue

        matched_mobile = match_atoms_by_proximity(ref_atoms, mobile_atoms)

        if len(matched_mobile) != len(ref_atoms):
            print("Could not match all atoms. Skipping pose.")
            continue

        sup = Superimposer()
        sup.set_atoms(ref_atoms, matched_mobile)
        rmsds.append(sup.rms)

    return rmsds

def dock_single_point(center):
    v = Vina(sf_name='vina')
    v.set_receptor('8IBT_receptor.pdbqt')
    v.set_ligand_from_file('LNT_ligand.pdbqt')

    # Use original clean PDB ligand structure as reference
    ref_structure = load_structure_from_file(REFERENCE_LIGAND_PATH)

    box_size = [18, 10, 13]
    v.compute_vina_maps(center=center, box_size=box_size)
    v.dock(exhaustiveness=32, n_poses=10)

    output_file = f'ligand_docked_{center[0]}_{center[1]}_{center[2]}.pdbqt'
    v.write_poses(output_file, n_poses=10, overwrite=True)

    docked_structures = []
    with open(output_file, 'r') as f:
        pose_lines = []
        for line in f:
            if line.startswith("MODEL"):
                pose_lines = []
            elif line.startswith("ENDMDL"):
                atom_lines = [l for l in pose_lines if l.startswith("ATOM") or l.startswith("HETATM")]
                structure = get_structure_from_atoms(atom_lines)
                docked_structures.append(structure)
            else:
                pose_lines.append(line)

    energies = v.energies()
    affinities = energies[:, 0].tolist() if energies.size > 0 else []

    rmsd_values = calculate_rmsd_with_superimposer(ref_structure, docked_structures)

    result = [
        datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(center),
        f"{np.mean(affinities):.2f}" if affinities else "No Pose",
        f"{np.std(affinities):.2f}" if affinities else "N/A",
        ', '.join(f"{a:.2f}" for a in affinities) if affinities else "N/A",
        ', '.join(f"{r:.2f}" for r in rmsd_values) if rmsd_values else "N/A"
    ]

    return result

grid_center = (75, -30, -60)
print(f"Starting docking for grid center: {grid_center}")
start_time = datetime.now()
result = dock_single_point(grid_center)
end_time = datetime.now()
print(f"Execution time: {end_time - start_time}")

with open('affinity_results_single.txt', 'w') as result_file:
    result_file.write('Timestamp\tGrid Center\tMean Affinity\tStd Dev\tAll Affinities\tRMSDs per Pose\n')
    result_file.write('---------\t-----------\t--------------\t--------\t----------------\t--------------\n')
    result_file.write('\t'.join(map(str, result)) + '\n')

df = pd.DataFrame([result], columns=["Timestamp", "Grid Center", "Mean Affinity", "Std Dev", "Affinities", "RMSDs"])
df.to_csv("affinity_results_single.csv", index=False)

print("Docking completed! Results saved in 'affinity_results_single.txt' and 'affinity_results_single.csv'")
