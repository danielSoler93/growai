import argparse
import growai.analysis.validity as vl
import growai.analysis.drug_like as dl
import growai.analysis.chemical_space as cs

def parse_args(parser):
    parser.add_argument('--sdf_files',  nargs="+", type=str, help='sdf files with the generated molecules', required=True)
    parser.add_argument('--druglikeness', action="store_true", help='Run druglikeness analysis')
    parser.add_argument('--validity', action="store_true", help="Run validity analysis")
    parser.add_argument('--pca', action="store_true", help="Retrieve the explored chemical space")



def main(sdf_files, druglikeness, validity, pca, xlabel, ylabel, xlim, ylim, title, output,
    growing_sites, iterations):
    if druglikeness:
        dl.plot_druglikeness(sdf_files, xlabel, ylabel, xlim, ylim, title, output)
    if validity:
        vl.get_validity(sdf_files, growing_sites, iterations)
    if pca:
        cs.plot_chemical_space(sdf_files, xlabel, ylabel,
            title, output, ylim, xlim)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='GrowingAI', conflict_handler='resolve')
    parse_args(parser)
    dl.parse_args(parser)
    vl.parse_args(parser)
    cs.parse_args(parser)
    args = parser.parse_args()
    main(args.sdf_files, args.druglikeness, args.validity, args.pca, args.xlabel, args.ylabel, args.xlim,
        args.ylim, args.title, args.output, args.growing_sites, args.iterations)

