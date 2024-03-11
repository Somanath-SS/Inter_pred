from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


def parse_ligand_pdb(file_path, ligand_chain='A'):
    ligand_atoms = []

    with open(file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("HETATM") and line[21] == ligand_chain:
                atom_type = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # Use the atom mapping feature to ensure correct atom types
                mapped_atom_type = f"{atom_type}_{len(ligand_atoms) + 1}"

                ligand_atoms.append((mapped_atom_type, (x, y, z)))

    return ligand_atoms


# Replace 'your_pdb_file.pdb' with the actual path to your PDB file
pdb_file_path = "E:/THESIS/BACK_END/4MS4-20240302T114112Z-001/4MS4/1.pdb"

# Parse ligand atoms with chain identifier 'A'
ligand_atoms = parse_ligand_pdb(pdb_file_path, ligand_chain='A')

# Create a molecule from ligand atoms
ligand_molecule = Chem.MolFromMolBlock('')
for atom_type, coords in ligand_atoms:
    atom = Chem.Atom(atom_type)
    ligand_molecule.AddAtom(atom)

# Add implicit hydrogens and generate a 2D structure diagram
ligand_molecule = Chem.AddHs(ligand_molecule)
AllChem.Compute2DCoords(ligand_molecule)

# Save the diagram as an image
img = Draw.MolToImage(ligand_molecule, size=(300, 300))
img.save("ligand_diagram.png")
img.show()
