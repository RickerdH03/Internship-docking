
# Built-in modules
import re
import subprocess

# Third-party modules
import numpy as np
from vina import Vina 

# TODO: score return several enervy values
# TODO: optimize poses

class AutoDockVina:

    def __init__(self, receptor_file: str, ligand_file: str):
        self.receptor_file = receptor_file
        self.ligand_file = ligand_file
        self.vina = Vina(sf_name='vina')
        self.vina.set_receptor(receptor_file)
        self.vina.set_ligand_from_file(ligand_file)

    def dock(
            self, 
            center: list, 
            box_size: list, 
            exhaustiveness: int = 8, 
            n_poses: int = 1,
            output_name: str = 'output',
            ):
        # Run docking
        self.vina.compute_vina_maps(center=center, box_size=box_size)
        print(f'Docking starting...')
        self.vina.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        
        # Manage results
        self.output_file = f'{output_name}.pdbqt'
        self.vina.write_poses(self.output_file, n_poses=n_poses, overwrite=True)
        print(f'Docking finished. Results saved to {self.output_file}.')
        #print(f'Spliting poses...')
        #cmd = f'obabel -ipdbqt {self.output_file}  -opdbqt -O {output_name}_*.pdbqt'
        #subprocess.run(cmd, shell=True)

        # Return energy values
        energies = self.vina.energies()
        return energies[:, 0]
    
    def report(self):
        
        # RMSD
        print(f'Calculating RMSD...')
        cmd = f'obrms -f {self.ligand_file} {self.output_file} -m'
        rmsds = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        rmsds = re.findall(r'([-+]?\d*\.\d+)[\n$]', rmsds.stdout)
        rmsds = np.array([float(rmsd) for rmsd in rmsds])

        # Energy
        print(f'Calculating energy...')
        energies = self.vina.energies()
        energies = energies[:, 0]

        return {
            'rmsd': rmsds,
            'energy': energies
        }


        
    
    

if __name__ == '__main__':
    receptor_file = '8IBT_receptor.pdbqt'
    ligand_file = 'LNT_ligand.pdbqt'
    
    vina = AutoDockVina(receptor_file, ligand_file)
    
    center = [75, -30, -60]
    box_size = [18, 10, 13]
    
    vina.dock(center=center, box_size=box_size, exhaustiveness=32, n_poses=10) 
    report = vina.report()
    print(report)
