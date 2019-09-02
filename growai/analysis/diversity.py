import os
from rdkit import Chem
from sklearn.decomposition import PCA
import visualization as vs
import molecule as ml


template_name = "round"
iterations = 10
initial_molecule = "STR_lig.sdf"

all_fingerprints = []
labels = []

sdf_file = initial_molecule
compound = next(Chem.SDMolSupplier(sdf_file))
m = ml.Molecule(compound)
all_fingerprints.append(m.get_morgan_fingerprint())
labels.append(True)

for i in range(iterations):
    print(template_name, i)
    name = template_name + str(i)
    sdf_file = os.path.join(name, name + ".sdf")
    suppl = Chem.SDMolSupplier(sdf_file)
    for molecule in suppl:
        if molecule:
            m = ml.Molecule(molecule)
            all_fingerprints.append(m.get_morgan_fingerprint())
            labels.append(False)
print(len(all_fingerprints))
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(all_fingerprints)
vs.plot(principalComponents[:,0], principalComponents[:, 1], labels, true_false=True)



