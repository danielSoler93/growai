import glob
from rdkit import Chem

FILES = "*/round*"


files = glob.glob(FILES)

all_valid_molecules = []
for file in files:
    molecules = Chem.SDMolSupplier(file)
    valid_molecules = [m for m in molecules if m]
    all_valid_molecules.extend(valid_molecules)

print("Molecules grown", len(all_valid_molecules))
