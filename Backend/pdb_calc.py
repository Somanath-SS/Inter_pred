import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import Bio.PDB


def parse_pdb(file_path):
    """
    Parse a PDB file and return the protein and ligand atoms.
    

    Args:
        file_path (str): The path to the PDB file to parse.

    Returns:
        _type_: _description_
    """
    protein_atoms = []
    ligand_atoms = []

    amino_acids = {}
    current_chain = None

    with open(file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                atom_type = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                chain = line[21]
                residue_number = int(line[22:26])
                residue_name = line[17:20]

                if chain != current_chain:
                    amino_acids = {}
                    current_chain = chain

                if chain not in amino_acids:
                    amino_acids[chain] = {}

                amino_acids[chain][residue_number] = residue_name

                protein_atoms.append(
                    (
                        atom_type,
                        x,
                        y,
                        z,
                        chain,
                        residue_number
                        )
                    )

            elif line.startswith("HETATM"):
                atom_type = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                ligand_atoms.append((atom_type, x, y, z))
    return protein_atoms, ligand_atoms, amino_acids

def calculate_distances(protein_atoms, ligand_atoms):
    """_summary_

    Args:
        protein_atoms (_type_): _description_
        ligand_atoms (_type_): _description_

    Returns:
        _type_: _description_
    """
    distances = []

    for protein_atom in protein_atoms:
       # print("Example protein atom:", protein_atom)  # Debug print
        for ligand_atom in ligand_atoms:
           # print("Example ligand atom:", ligand_atom)  # Debug print
            dist = np.sqrt((protein_atom[1] - ligand_atom[1]) ** 2 +
                           (protein_atom[2] - ligand_atom[2]) ** 2 +
                           (protein_atom[3] - ligand_atom[3]) ** 2)
            distances.append(
                    (
                    protein_atom[0],
                    ligand_atom[0],
                    dist,
                    protein_atom[4],
                    protein_atom[5]
                    )
                )

    return distances

def sort_distances(distances):
    """_summary_

    Args:
        distances (_type_): _description_

    Returns:
        _type_: _description_
    """
    return sorted(distances, key=lambda x: x[2])

def find_hydrogen_bonds(sorted_distances):
    """_summary_

    Args:
        sorted_distances (_type_): _description_

    Returns:
        _type_: _description_
    """
    hydrogen_bonds = []

    for protein_atom_type, ligand_atom_type, dist, chain, residue_number \
        in sorted_distances:
        if (protein_atom_type[0] in ['N', 'O', 'S']  \
            and ligand_atom_type[0] in ['N', 'O', 'S']) or \
                (ligand_atom_type[0] in ['N', 'O', 'S'] \
                    and protein_atom_type[0] in ['N', 'O', 'S']):
            if 2.7 <= dist <= 3.3:  # Adjust the threshold distance 
                                    # for hydrogen bonds
                hydrogen_bonds.append(
                    (chain,
                     residue_number,
                     protein_atom_type,
                     ligand_atom_type,
                     dist
                     )
                    )
                
    if hydrogen_bonds:
        print("\nPossible Hydrogen Bond Interactions:")
        for bond in hydrogen_bonds:
            print(bond)
    else:
        print("\nNo possible Hydrogen Bond Interactions found.")

    return hydrogen_bonds

def visualize_ligand_sticks_with_labels(
    ligand_atoms,
    hydrogen_bonds,
    protein_atoms,
    amino_acids
    ):
    """_summary_

    Args:
        ligand_atoms (_type_): _description_
        hydrogen_bonds (_type_): _description_
        protein_atoms (_type_): _description_
        amino_acids (_type_): _description_
    Returns:
        _type_: _description_
    """
    # Create a Py3Dmol view
    view = py3Dmol.view(width=800, height=600)

    # Create a PDB file string for ligand atoms
    ligand_pdb_string = ""
    for i, (atom_type, x, y, z) in enumerate(ligand_atoms):
        ligand_pdb_string += (
        f"ATOM  {i+1:<5} {atom_type:<4} LIG     1    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}\n"
                                )
    # Add ligand atoms to the view
    view.addModel(ligand_pdb_string, "pdb")

    # Set style to stick representation for ligand atoms
    view.setStyle({'model': 0}, {'stick': {}})

    # Add labels to ligand atoms
    for i, (atom_type, x, y, z) in enumerate(ligand_atoms):
        view.addLabel(atom_type, {'position': {'x': x, 'y': y, 'z': z}})

    # Create a PDB file string for protein atoms
    protein_pdb_string = ""
    for i, (atom_type, x, y, z, chain, residue_number) \
        in enumerate(protein_atoms):
        protein_pdb_string += (
            f"ATOM  {i+1:<5} {atom_type:<4} PRO     1"
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}\n"
                               )

    # Add protein atoms to the view
    view.addModel(protein_pdb_string, "pdb")

    # Set style to sphere representation for protein atoms
    view.setStyle({'model': 1}, {'sphere': {'color': 'red'}})

    # Add labels to protein atoms
    for chain, residue_number, _, _, _ in hydrogen_bonds:
        amino_acid = amino_acids.get(chain, {}).get(
            residue_number, "Unknown"
            )
        # Find the coordinates of the protein atom involved in the hydrogen bond
        protein_atom_coords = next(
            (x, y, z) for _, x, y, z, ch, rn  \
                in protein_atoms if ch == chain and rn == residue_number
                )
        # Label the protein atom with amino acid and residue number
        label_text = f"{amino_acid} {residue_number}({chain})"
        view.addLabel(
            label_text, {
                          'position': 
                              {
                              'x': protein_atom_coords[0],
                              'y': protein_atom_coords[1],
                              'z': protein_atom_coords[2]
                              }
                          }
                      )

        # Find the ligand atom coordinates involved in the hydrogen bond
        ligand_atom_coords = None
        for _, x, y, z in ligand_atoms:
            if _ == "N1" or _ == "O2":
                ligand_atom_coords = (x, y, z)
                break

        # Connect protein atom to ligand with yellow dotted line
        view.addLine(
                        {
                        'start': {
                            'x': protein_atom_coords[0],
                            'y': protein_atom_coords[1],
                            'z': protein_atom_coords[2]
                                },
                        'end': {
                            'x': ligand_atom_coords[0],
                            'y': ligand_atom_coords[1],
                            'z': ligand_atom_coords[2]
                                },
                        'color': 'yellow', 
                        'dashed': True
                        }
                    )

    # Set background color to black
    view.setBackgroundColor('black')

    # Zoom to fit the structure
    view.zoomTo()

    # Show the 3D view
    view.show()


if __name__ == "__main__":
    protein_atoms, ligand_atoms, amino_acids = parse_pdb(r"..\sample_data\first.pdb")
    distances = calculate_distances(protein_atoms, ligand_atoms)
    sorted_distances = sort_distances(distances)
    hydrogen_bonds = find_hydrogen_bonds(sorted_distances)
    visualize_ligand_sticks_with_labels(ligand_atoms, hydrogen_bonds, protein_atoms, amino_acids)
    # Finding and printing possible hydrogen bonds
    hydrogen_bonds = find_hydrogen_bonds(sorted_distances)
    if hydrogen_bonds:
        print("\nPossible Hydrogen Bond Interactions:")
        for chain, residue_number, protein_atom_type, ligand_atom_type, \
            dist in hydrogen_bonds:
            amino_acid = amino_acids.get(
                chain, {}).get(residue_number, "Unknown")
            print(f"{amino_acid} {residue_number}({chain})"
                  f"{protein_atom_type}(donor) ------ {ligand_atom_type}"
                  f"(acceptor)   {dist}")
    else:
        print("\nNo possible Hydrogen Bond Interactions found.")
        
    visualize_ligand_sticks_with_labels(
        ligand_atoms,
        hydrogen_bonds,
        protein_atoms,
        amino_acids
        )

