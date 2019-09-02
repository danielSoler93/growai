from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import qed


class Molecule():


    def __init__(self, molecule):
        self._mol = molecule
    
    def get_morgan_fingerprint(self):
        self._morgan_fingerprint = AllChem.GetMorganFingerprintAsBitVect(self._mol, 2)
        self.morgan_fingerprint = [bit for bit in self._morgan_fingerprint.ToBitString()]
        return self.morgan_fingerprint, self._morgan_fingerprint


    def get_daylight_fingerprint(self):
        self._daylight = FingerprintMols.FingerprintMol(self._mol)
        self.daylight = [bit for bit in self._daylight.ToBitString()]
        return self.daylight, self._daylight

    def get_tanimoto(self, mol_ref):
        _, daylight_fp1 = self.get_daylight_fingerprint()
        _, daylight_fp2 = self.get_daylight_fingerprint()
        self.tanimoto = DataStructs.FingerprintSimilarity(daylight_fp1, daylight_fp2) 
        return self.tanimoto

    def get_qed(self):
        self.qed = qed(self._mol) 
        return self.qed
