import argparse
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import growai.analysis.helpers as hp
import numpy as np


def plot_chemical_space(files, xlabel="First component PCA", ylabel="Second component PCA", 
    title="PCA", output="pca.png", ylim=None, xlim=None):
    all_molecules = hp.retrieve_mols_from_files(files)
    fingerprints = [m.get_morgan_fingerprint()[1] for m in all_molecules if m]
    pca = PCA(n_components=2)
    pca = pca.fit_transform(fingerprints)
    plt.scatter(pca[1:,0], pca[1:,1], c="green")
    plt.scatter(pca[0,0], pca[0,1], c="black")
    plt.xlabel(xlabel)
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    plt.ylabel(ylabel)
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
    plot_chemical_space(args.sdf_files, args.xlabel, args.ylabel, args.title, args.output, args.ylim, args.xlim)

