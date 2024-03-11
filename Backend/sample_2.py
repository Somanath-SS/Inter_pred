import numpy as np
import pandas as pd


def parse_pdb(file_path):
    protein_atoms = []
    ligand_atoms = []

    with open(file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                atom_type = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                protein_atoms.append((atom_type, x, y, z))

            elif line.startswith("HETATM"):
                atom_type = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                ligand_atoms.append((atom_type, x, y, z))

    return protein_atoms, ligand_atoms


def calculate_distances(protein_atoms, ligand_atoms):
    distances = []

    for protein_atom in protein_atoms:
        for ligand_atom in ligand_atoms:
            dist = np.sqrt((protein_atom[1] - ligand_atom[1]) ** 2 +
                           (protein_atom[2] - ligand_atom[2]) ** 2 +
                           (protein_atom[3] - ligand_atom[3]) ** 2)
            distances.append((protein_atom[0], ligand_atom[0], dist))

    return distances


def sort_distances(distances):
    return sorted(distances, key=lambda x: x[2])


def write_to_excel(sorted_distances, output_file='distances.xlsx'):
    df = pd.DataFrame(sorted_distances, columns=['Protein Atom', 'Ligand Atom', 'Distance'])
    df.to_excel(output_file, index=False)
    print(f"Distances written to {output_file}")


def find_hydrogen_bonds(sorted_distances, threshold_distance=2.5, num_bonds_to_check=20):
    hydrogen_bonds = []

    for i in range(min(num_bonds_to_check, len(sorted_distances))):
        protein_atom_type, ligand_atom_type, dist = sorted_distances[i]

        # Check for possible hydrogen bonds, considering pairs like "H O"
        if (protein_atom_type[0] == 'H' and ligand_atom_type[0] in ['O', 'N', 'S']) or \
                (ligand_atom_type[0] == 'H' and protein_atom_type[0] in ['O', 'N', 'S']):
            hydrogen_bonds.append((protein_atom_type, ligand_atom_type, dist))

    return hydrogen_bonds



pdb_file_path = "E:/THESIS/BACK_END/4MS4-20240302T114112Z-001/4MS4/1.pdb"

protein_atoms, ligand_atoms = parse_pdb(pdb_file_path)
distances = calculate_distances(protein_atoms, ligand_atoms)
sorted_distances = sort_distances(distances)

#sorted distances to an Excel file
write_to_excel(sorted_distances, output_file='sorted_distances.xlsx')

# Finding and printing possible hydrogen bonds
hydrogen_bonds = find_hydrogen_bonds(sorted_distances)
if hydrogen_bonds:
    print("\nPossible Hydrogen Bond Interactions:")
    for protein_atom_type, ligand_atom_type, dist in hydrogen_bonds:
        print(f"{protein_atom_type}(donor) ------ {ligand_atom_type}(acceptor) \t {dist}(distance)")
else:
    print("\nNo possible Hydrogen Bond Interactions found.")
