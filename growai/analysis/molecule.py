from rdkit.Chem import AllChem


class Molecule():


    def __init__(self, molecule):
        self._mol = molecule
    
    def get_morgan_fingerprint(self):
        fp_rdkit = AllChem.GetMorganFingerprintAsBitVect(self._mol, 2).ToBitString()
        fp = [bit for bit in fp_rdkit]
        self.morgan_fingerprint = fp
        return fp
