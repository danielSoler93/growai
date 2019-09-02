from rdkit import Chem

FILE = "fragments.sdf"


suppl = Chem.SDMolSupplier(FILE)

molecules = []
for m in suppl:
    ring_info = m.GetRingInfo()
    try:
        natoms_in_ring = 0
        for number_of_ring_atoms in ring_info.AtomRings():
            if number_of_ring_atoms > natoms_in_ring:
                natoms_in_ring = number_of_ring_atoms
    except IndexError:
        natoms_in_ring = 0
    if natoms_in_ring <= 6:
        molecules.append(m)
print("Synthtically feasible molecules", len(molecules))
