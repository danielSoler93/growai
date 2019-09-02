from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from rdkit.Chem import AllChem
from tqdm import tqdm
import growai.analysis.molecule as ml




def retrieve_mols_from_files(files):

    print("Retrieving molecules")
    all_molecules = []
    for file in tqdm(files):
        molecules = Chem.SDMolSupplier(file)
        for molecule in molecules:
            if molecule:
                mol_python = ml.Molecule(molecule)
                all_molecules.append(mol_python)
    return all_molecules


