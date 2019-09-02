import argparse
from rdkit import Chem
import glob
import matplotlib.pyplot as plt
from rdkit.Chem.Descriptors import qed
from rdkit.Chem import rdBase, RDConfig
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if
DrawingOptions.bondLineWidth=1.8
import growai.analysis.helpers as hp



def plot_druglikeness(files, xlabel="Algorithm", ylabel="QED", xlim=None, ylim=None,
    title="Drug likeness", output="druglikeness.png"):

    all_molecules = hp.retrieve_mols_from_files(files)
    
    QED = []
    for m in all_molecules:
      try:
        QED.append(m.get_qed())
      except:
        pass
    
    
    #plt.hist(QED)
    plt.boxplot(QED)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if ylim:
        plt.ylim(ylim)
    if xlim:
        plt.xlim(xlim)
    plt.title(title)
    plt.savefig(output)


def parse_args(parser):
    parser.add_argument('--sdf_files',  nargs="+", type=str, help='sdf files with the generated molecules', required=True)
    parser.add_argument('--xlabel', type=str, help='label of the x axis of the plot', default="QED")
    parser.add_argument('--ylabel', type=str, help='label of the y axis of the plot', default="Algorithm")
    parser.add_argument('--ylim', nargs="+", type=float, help='Limit of the y axis. i.e 0 1', default=None)
    parser.add_argument('--xlim', nargs="+", type=float, help='Limit of the x axis. i.e 0 1', default=None)
    parser.add_argument('--title', type=str, help='title of the plot', default="drug likeness")
    parser.add_argument('--output', type=str, help='output file name (png)', default="druglikeness.png")
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='GrowingAI')
    parse_args(parser)
    args = parser.parse_args()
    plot_druglikeness(args.sdf_files, args.xlabel, args.ylabel, args.xlim, args.ylim, args.title, args.output)

