from rdkit import Chem
from rdkit.Chem import Draw

def extract_ligand_smiles_from_pdb(pdb_file_path):
    ligand_block_started = False
    ligand_atoms = []

    with open(pdb_file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('HETATM') and 'UNK' in line:
                ligand_block_started = True
                # Extract atom type from the end of the line
                atom_type = line.strip().split()[-1]
                ligand_atoms.append(atom_type)
            elif ligand_block_started and line.startswith('TER'):
                break

    ligand_smiles = ''.join(ligand_atoms)
    return ligand_smiles

def draw_molecule_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        img.show()

# Replace 'your_pdb_file.pdb' with the actual path to your PDB file
pdb_file_path = "E:/THESIS/BACK_END/4MS4-20240302T114112Z-001/4MS4/1.pdb"

ligand_smiles = extract_ligand_smiles_from_pdb(pdb_file_path)
print('Ligand SMILES:', ligand_smiles)

draw_molecule_from_smiles(ligand_smiles)
