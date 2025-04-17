# Define docking function
def dock_single_point(center):
    v = Vina(sf_name='vina')
    v.set_receptor('8IBT_receptor.pdbqt')
    v.set_ligand_from_file('LNT_ligand.pdbqt')

    box_size = [18, 10, 13]
    affinities = []
    rmsd_values = []

    for repeat in range(1):  # Run 1 docking repeat
        v.compute_vina_maps(center=center, box_size=box_size)
        v.dock(exhaustiveness=32, n_poses=10)

        combined_output_file = f'ligand_docked_{center[0]}_{center[1]}_{center[2]}_rep{repeat + 1}.pdbqt'
        v.write_poses(combined_output_file, n_poses=10, overwrite=True)

        # Split individual poses to separate files
        with open(combined_output_file, 'r') as f:
            pose_idx = 1
            pose_lines = []
            for line in f:
                if line.startswith("MODEL"):
                    pose_lines = [line]
                elif line.startswith("ENDMDL"):
                    pose_lines.append(line)
                    pose_file = f'ligand_docked_{center[0]}_{center[1]}_{center[2]}_pose{pose_idx}.pdbqt'
                    with open(pose_file, 'w') as pf:
                        pf.writelines(pose_lines)
                    pose_idx += 1
                else:
                    pose_lines.append(line)

        energy = v.energies()
        if energy.size > 0:
            affinities.append(energy[0, 0])
