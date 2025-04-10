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
### RMSDtest.py
Simplified version of parralelnewgrid.py with an RMSD calculation added. It runs only 1 grid center (middle of active site). The RMSD is calculated using Biopython. Affinity values are not properly put into the output txt file but doesn't matter for now.
### pose.py
Short function for some docking settings and splitting the poses into seperate files. Script uses combined_output_file.Adjust rest of script accordingly.

## Expected output
The expected out is an RMSD between 2 and 12 across the different poses. RMSD is not dependent on the ranking of the poses- pose 8 can have RMSD 2 while pose 1 can have RMSD 10. If the RMSD's are all below 2 and/or very close to eachother (+/- 1) something has gone wrong. If the RMSD is around 100 you must first align the poses to the ligand file.

