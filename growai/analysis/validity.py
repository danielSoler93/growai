import argparse
import glob
from rdkit import Chem




def get_validity(files, growing_sites, iterations):

    all_valid_molecules = []
    for file in files:
        molecules = Chem.SDMolSupplier(file)
        valid_molecules = [m for m in molecules if m]
        all_valid_molecules.extend(valid_molecules)
    
    print("Theoretical number of molecules", 4**iterations*growing_sites)
    print("Valid molecules", len(all_valid_molecules))
    print("Sucess rate", len(all_valid_molecules)/4**iterations*100/growing_sites, "%")



def parse_args(parser):
    parser.add_argument('--sdf_files',  nargs="+", type=str, help='sdf files with the generated molecules', required=True)
    parser.add_argument('--growing_sites', type=int, help='Growing sites identified by the software in the original molecule', default=1)
    parser.add_argument('--iterations', type=int, help='Growing iterations performed', default=5)
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='GrowingAI')
    parse_args(parser)
    args = parser.parse_args()
    get_validity(args.sdf_files, args.growing_sites, args.iterations)

