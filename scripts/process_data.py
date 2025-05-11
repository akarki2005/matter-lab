# Script to calculate polarizability tensor of molecules given SMILES structure as input
# Written by Aayush Karki

import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem

# Read input xlsx file into a pandas dataframe
df = pd.read_excel('../database/pre-screening-data.xlsx')

# Print head (testing, remove later)
# print(df.head())

# Process data

# TODO

# Algorithm for converting a SMILES string to an xyz file
def smiles_to_xyz(smiles: str, output_file: str) -> None:
    """
    Parameters
    ----------
    smiles : str
        The SMILES representation of the molecule
    output_file : str
        The path to the file that you want to store the xyz of the molecule in

    Returns 
    -------
    Nothing; mutates the output file to contain the xyz representation of the molecule

    Raises
    ------
    FileNotFoundError if the path to the file is invalid
    ValueError if smiles input is invalid, or some rdkit operation fails during execution

    """
    
    # convert smiles to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES string: {smiles}")

    # add hydrogens
    mol = Chem.AddHs(mol)

    # construct 3D structure for the molecule
    result = AllChem.EmbedMolecule(mol)
    if result == -1: # error
        raise ValueError("failed to generate 3D structure")

    # optimize spatial geometry
    result = AllChem.UFFOptimizeMolecule(mol)
    if result == -1:
        raise ValueError("failed to optimize geometry")

    conformer = mol.GetConformer()

    with open(output_file, "w") as file:
        # write in initial information
        file.write(f"{mol.GetNumAtoms()} \n") # write number of atoms
        file.write("Molecule Name (placeholder) \n") # write name of molecule (not complete, placeholder name for now)
        # write in positions of all atoms
        for atom in mol.GetAtoms():
            # get atom position (x, y, z coords)
            index = atom.GetIdx()
            position = conformer.GetAtomPosition(index)
            # write atom position to file, to 5 decimal places
            file.write(f"{atom.GetSymbol()} {position.x:.5f} {position.y:.5f} {position.z:.5f} \n")

smiles_to_xyz('Ic1ccc(I)cc1', '../database/test_output.txt')
