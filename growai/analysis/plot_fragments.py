from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if
DrawingOptions.bondLineWidth=1.8


FILE = "fragments.sdf"

suppl = [m for m in Chem.SDMolSupplier(FILE)]


for m in suppl: tmp=AllChem.Compute2DCoords(m)


img = Draw.MolsToGridImage([m for m in suppl], molsPerRow=40)

img.save("image.png")

