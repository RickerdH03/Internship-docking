# Internship-docking

## File description
### 8IBT_receptor.pdbqt
Receptor input for vina. It is the B. infantis beta-galactosidase as a monomer with the complexed LNT removed and converted to .pdbqt
### LNT_ligand.pdbqt
Ligand input for vina. The complexed LNT was extracted into a seperate file and converted from .pdb to .pdbqt
### 20250213-LNT from PDB.pdb
Complexed LNT extracted into seperate file
### parralelnewgrid.py
Parralel running of docking across several grid centers within the active site of the B.infantis beta-galactosidase (8IBT)
### RMSDalright.py
Simplified version of parralelnewgrid.py with an RMSD calculation added. It runs only 1 grid center (middle of active site). The RMSD is calculated using OpenBabel. Affinity values are not properly put into the output txt file but doesn't matter for now.
### RMSDtest
Simplified version of parralelnewgrid.py with an RMSD calculation added. It runs only 1 grid center (middle of active site). The RMSD is calculated using Biopython. Affinity values are not properly put into the output txt file but doesn't matter for now.
