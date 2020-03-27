import os
from rdkit import Chem
from rdkit import rdBase
from tqdm import tqdm
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if
import scipy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
from rdkit.Chem import rdMolAlign
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, median_absolute_error
from sklearn.model_selection import cross_val_predict, cross_val_score
from sklearn.metrics import confusion_matrix
import tensorflow as tf
from keras.optimizers import RMSprop
from keras.models import Sequential
from keras.layers import Dense, Activation
from keras.layers import LSTM, Input, Flatten
from keras.layers import Dropout,  BatchNormalization
from keras.layers import Embedding, Bidirectional
from keras.utils import to_categorical
from keras.preprocessing.sequence import pad_sequences
from tqdm import tqdm
import subprocess
import matplotlib.pyplot as plt
import collections
import pandas as pds
import prody as pd
import time
import glob
import argparse
import numpy as np
import growai.model.docking.glide as gl
import growai.model.helpers.helpers as hp
import growai.analysis as an

# Deactive tnesorflow warnings
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

# Deactivate prody warnings
pd.confProDy(verbosity="none")

# Disable Rdkit stdout
rdBase.DisableLog('rdApp.error')

#Harcode RDKIT Data
#RDKIT=os.path.join(cs.SCHRODINGER, "mmshare-v*/data/RDKit/Data/")
RDKIT=""

#Print out Error if RDKIT not recognize
try:
    fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
except IOError:
    print("If error: not RDDataDir specified on installation. Manually change the RDKIT variable on this file {}".format(os.__file__))
    fdefName = os.path.join(RDKIT,'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    

# Classifiersi:

###FIRST MODEL TO STACK####
from sklearn.svm import SVR
CLF1 = SVR(gamma=0.001, C=10.0, kernel="rbf")
param_grid = [ { 'C': [1, 10, 100, 1000], "kernel": ["linear"]},
               { 'C': [1, 10, 100, 1000], "kernel": ["rbf"], "gamma": [0.001, 0.0001]},]
####SECOND MODEL TO STACK####
from sklearn import linear_model
CLF2 = linear_model.Ridge(alpha=100, normalize=False, fit_intercept=True)
param_grid = {"alpha":[1000, 500, 100, 25, 10, 5, 1, 0.8, 0.5, 0.3, 0.1 , 0.01]}

####THIRD MODEL TO STACK###
from sklearn.neighbors import KNeighborsRegressor
CLF3 = KNeighborsRegressor(n_neighbors=3)

####FORTH MODEL TO STACK###
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
kernel = DotProduct() + WhiteKernel()
CLF4 = GaussianProcessRegressor(kernel=kernel, random_state=0)

###STACK MODEL####
CLF5 = SVR(gamma=0.001, C=10.0, kernel="rbf")

###ASSESMENT MODEL###
CLF6 = SVR(gamma=0.001, C=10.0, kernel="rbf")

# Growing Vocabulary
VOC = {'!': 0,
 '(': 3,
 ')': 5,
 '1': 8,
 '2': 11,
 '3': 12,
 '=': 4,
 'C': 1,
 'N': 6,
 'O': 2,
 'c': 10,
 'n': 9,
 '|': 7}

# Aminoacid Vocabulary
VOC_AMINOACIDS = ["ARG", "LYS", "GLN", "ASN", "HIS", "TYR", "TRP", "SER", "THR", "PHE", "PRO", "ILE", "VAL", "LEU", \
      "MET", "CYS",  "ALA", "GLY", "ASP", "GLU" ]

# Colum names
MACCS_COLUMNS = ["MACS_{}".format(i) for i in range(1, 167)]
COLUMNS= [ "n_carbon", "n_nitrogen", "n_oxigen", 'n_Donor','n_Acceptor', 'Aromatic', 'Hydrophobe', \
      'n_atoms'] + MACCS_COLUMNS + [ "rotatable_bonds", "ring", "HD", "HA", "n_het", "n_amide", "sp3", "MolWT", "x", "y", "z", 'closest_residue', 'average_residue_dist', 'num_residues_at_5A',\
      "ARG", "LYS", "GLN", "ASN", "HIS", "TYR", "TRP", "SER", "THR", "PHE", "PRO", "ILE", "VAL", "LEU", \
      "MET", "CYS",  "ALA", "GLY",  "ASP", "GLU" ]

def split_complex(pdb, resname, output_ligand=None, output_rec=None):
    with open(pdb, "r") as f:
        lines = f.readlines()
        ligand = [line for line in lines if line[17:20].strip() == resname ]
        receptor = [line for line in lines  if (line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER")) and line[17:20].strip() != resname ]
    if not output_ligand:
        output_ligand = resname + "_lig.pdb"
    if not output_rec:
        output_rec = resname + "_rec.pdb"
    with open(output_ligand, "w") as fout:
        fout.write("".join(ligand))
    with open(output_rec, "w") as fout:
        fout.write("".join(receptor))
    return os.path.abspath(output_ligand), os.path.abspath(output_rec)

def pdb_to_sdf(ligand_pdb, schr="", output=None):
    mol = Chem.rdmolfiles.MolFromPDBFile(ligand_pdb, removeHs=False)
    if not output:
        output = os.path.splitext(ligand_pdb)[0] + ".sdf" 
    w = Chem.SDWriter(output)
    w.write(mol)
    #babel_loc = os.path.join(schr, "utilities/babel")
    #command = "{} -ipdb {} -osdf {}".format(babel_loc, ligand_pdb,  output)
    #subprocess.call(command.split())
    return os.path.abspath(output)

def grow_protocol(pdb, sdf, vocabulary, resname, max_len=20, iterations=10, grow=True, rank=True):
    print("\nBuild growing model")
    print("===================\n")
    model = create_model()
    
    print("\nBuild initial molec")
    print("===================\n")
    core_molecule = Chem.AddHs(next(Chem.SDMolSupplier(sdf)))
    core_smile = Chem.MolToSmiles(next(Chem.SDMolSupplier(sdf)))
    core_smile_stripped = core_smile.strip("[HH].").strip("([H])")
    core_smile_pattern, n_growing_positions = identify_spots_to_grow(core_smile_stripped)

    #Add bonded atoms to COLUMN
    stack_models = None
    models = []
    cold_start = True
    rsquares = []
    scores = []
    steps = []
    total_fragments = []
    number_of_total_fragments = 0
    for i in range(n_growing_positions):
        COLUMNS.insert(0, "bonded_to_atom{}".format(n_growing_positions-i))

    print("\nGrowing")
    print("=========\n")
    smiles_current_generation = "!"
    n_mol = 0
    final_smiles = []
    features = pds.DataFrame(columns=COLUMNS) 
    labels = []
    use_docking_model = False
    previous_score = 100000

    if grow:

        for i in tqdm(range(iterations)):
          smiles_next_generation = []
          for j, smile in enumerate(smiles_current_generation):
            if len(smile) > max_len:
              continue
            else:
              predicted_smiles = grow_molecule(smile, vocabulary, model)
              for predicted_smile in predicted_smiles:
                predicted_smile = predicted_smile.strip("|")
                smiles_next_generation.append(predicted_smile)

          output_folder = 'round{}'.format(i)
          output_sdf = 'round{}.sdf'.format(i)

          print("\nOutput molecules to {}".format(os.path.join(output_folder, output_sdf)))
          print("========================\n")
          fragments = list(set([smile for smile in smiles_next_generation if Chem.MolFromSmiles(smile.strip("!"))]))
          total_fragments.extend(fragments)
          number_of_total_fragments += len(fragments)
          final_molecules, fragments, bonded_atom = combine_smiles(core_smile_pattern, fragments, n_growing_positions, n_atoms=i+1)
          if not os.path.exists(output_folder):
              os.mkdir(output_folder)
          with(hp.cd(output_folder)):
              w = Chem.SDWriter(output_sdf)
              molecules_threeD = []
              for m, f, b in zip(final_molecules, fragments, bonded_atom):
                  try:
                      m_h = Chem.AddHs(m)
                  except ValueError:
                      print("Skipped")
                      continue 
                  AllChem.Compute2DCoords(m_h)
                  m_h.SetProp("_Name","molecule{}".format(n_mol))
                  m_h.SetProp("fragment", f)
                  m_h.SetProp("bonded_atom", str(b))
                  n_mol += 1
                  w.write(m_h)
                  molecules_threeD.append(m_h)
          final_smiles.extend(smiles_current_generation)
          smiles_current_generation = smiles_next_generation
          print("{} fragments generated".format(number_of_total_fragments))
          print("{} molecules generated".format(len(set(final_smiles))))

    # Write sdf with fragments
    total_fragments_mols = [Chem.MolFromSmiles(smile.strip("!")) for smile in total_fragments]
    w = Chem.SDWriter("fragments.sdf")
    for m in total_fragments_mols:
        try:
            m_h = Chem.AddHs(m)
        except ValueError: 
            continue
        AllChem.Compute2DCoords(m_h)
        w.write(m_h)

    j = 0
    molecules = []
    folders = glob.glob("round*")
    #folders = glob.glob("round4")

    if rank:

        print("\tExtract molecules of {} files".format(len(folders)))
        for folder in folders:
            sdf_file = os.path.join(folder, folder + ".sdf")
            molecules.extend([m for m in tqdm(Chem.SDMolSupplier(sdf_file))])
        fp_ref = FingerprintMols.FingerprintMol(molecules[0])

        print("\nCompute Similarity")
        print("===================\n")
        similarity = np.array([ DataStructs.FingerprintSimilarity(fp_ref, FingerprintMols.FingerprintMol(m)) if m else 10 for m in tqdm(molecules)])
        molecules_idxs_sort_by_similarity = np.argsort(similarity).tolist()
        similarity = similarity.tolist()
        mol_ref = similarity[molecules_idxs_sort_by_similarity[0]]

        while True:

            # Create Folder
            if j > 100:
                sys.exit(-1)
            folder_name = "learning_{}".format(j)
            if not os.path.exists(folder_name):
                os.mkdir(folder_name)
            
            #Move inside Folder
            with(hp.cd(folder_name)):
                
                # If Predicting docking score (once learned)
                if use_docking_model:
                    features, features_test = retrieve_features(features, pdb, molecules, fragments, core_molecule, 
                        resname,  n_growing_positions, bonded_atom, labels=None)
                    final_features = features_test.copy()
                    for i, models_round in enumerate(models):
                        if i == len(models)-1:
                            intermidiate_features = features_test.copy()
                            model1, model2, model3, model4, model_final = models_round 
                            pred1 = model1.predict(scaler.fit_transform(features_test))
                            pred2 = model2.predict(scaler.fit_transform(features_test))
                            pred3 = model3.predict(scaler.fit_transform(features_test))
                            pred4 = model4.predict(scaler.fit_transform(features_test))
                            intermidiate_features["SVM".format(i)] = pred1
                            intermidiate_features["linear".format(i)] = pred2
                            intermidiate_features["KN".format(i)] = pred3
                            intermidiate_features["GP".format(i)] = pred4
                            final_features["final_model_{}".format(i)] = model_final.predict(scaler.fit_transform(intermidiate_features))
                        else:
                            final_features["final_model_{}".format(i)] = [0] * features_test.shape[0]
                    docking_score = model_to_docking.predict(scaler.fit_transform(final_features))
                    print(docking_score)
                    if os.path.exists(data_file):
                        # Score
                        data_file = an.analyze(glob.glob(output_docking_mae), filter=["Score"])
                        test_labels = retrieve_labels(data_file)
                        value = median_absolute_error(test_labels, docking_score)
                        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(test_labels, docking_score)
                        print("Score:", value)
                        print("r2", r_value)
                        scores.append(value)
                        rsquares.append(r_value)
                        steps.append(j)
                        # Plot
                        fig, ax = plt.subplots()
                        ax.plot(steps, rsquares)
                        ax.plot(steps, scores)
                        fig.savefig("rsquare_scores.png")
                        # Plot
                        fig, ax = plt.subplots()
                        ax.scatter(test_labels, docking_score)
                        fig.savefig("result.png")
                   
                # If learning
                else:
                    if cold_start:
                        # Elect structs to dock
                        number_of_structs = 1000
                        #number_of_structs = 20
                        indxs_to_pick = np.round(np.linspace(0, len(molecules_idxs_sort_by_similarity) - 1, number_of_structs)).astype(int)
                        indxs_to_pick = sorted(indxs_to_pick, reverse=True)
                        indxs_mol = [molecules_idxs_sort_by_similarity[idx] for idx in indxs_to_pick]

                        molecules_to_dock = []
                        for i, (idx, idx_remove) in enumerate(zip(indxs_mol, indxs_to_pick)):
                            if similarity[idx] <= 1:
                                #Save molec
                                m = molecules[idx]
                                m.SetProp("Similarity", str(similarity[idx]))
                                molecules_to_dock.append(molecules[idx])
                                #Remove molecules from eligible indexes
                                molecules_idxs_sort_by_similarity.pop(idx_remove)
                            else:
                                #If molec is wrong remove it from eligible
                                molecules_idxs_sort_by_similarity.pop(idx_remove)
                        cold_start = False
                    else:
                        new_molecules = []
                        for k, (m_ref, result) in tqdm(enumerate(zip(molecules_to_dock, results))):
                            limit_similar_mols = 15
                            number_of_similar_mols = 0
                            sim_ref = float(m_ref.GetProp("Similarity"))
                            if result > 1.5: 
                                #If the result differs in a docking unit choose similar molecules
                                idxs = np.argsort(np.abs(np.array(similarity) - sim_ref))
                                for i in idxs:
                                    if number_of_similar_mols > limit_similar_mols:
                                        break 
                                    if i in molecules_idxs_sort_by_similarity and molecules[i] and molecules[i] not in new_molecules:
                                        chosen_mol = molecules[i]
                                        chosen_mol.SetProp("Similarity", str(similarity[i]))
                                        new_molecules.append(chosen_mol)
                                        molecules_idxs_sort_by_similarity.remove(i)
                                        number_of_similar_mols += 1
                        molecules_to_dock = new_molecules
                        print("Number of Chosen Molecules", len(molecules_to_dock))
                                   

                    #Write input for docking
                    sdf = "structs_to_dock_{}.sdf".format(i)                
                    w = Chem.SDWriter(sdf)
                    w.SetProps([""])
                    for m in molecules_to_dock: w.write(m)
                    
                    #Dock
                    print("\nDocking")
                    print("===========\n")
                    output_docking_mae = "input__*__dock_lib.maegz" 
                    output_docking = "input__*__dock_subjob_poses.zip"
                    if not glob.glob(output_docking):
                        dock(pdb, sdf)
                    while not glob.glob(output_docking):
                        time.sleep(60)

                    print("\nFeature Extraction")
                    print("=====================\n")
                    data_file = an.analyze(glob.glob(output_docking_mae), filter=["Score"]) 
                    molecules_to_dock, test_labels = retrieve_labels(molecules_to_dock, data_file)
                    #test_labels = [-0, -1, -2, -3, -4, -0, -1, -2, -3, -4,-0, -1, -2, -3, -4,-0, -1, -2, -3, -4]
                    # Retrieve other info from final moleculs
                    fragments = [ m.GetProp("fragment") for m in molecules_to_dock ]
                    bonded_atom = [ m.GetProp("bonded_atom") for m in molecules_to_dock ]

                    #Build features
                    features, features_test, test_labels = retrieve_features(features, pdb, molecules_to_dock, fragments, core_molecule,
                        resname, n_growing_positions, bonded_atom, core_smile_pattern, labels=test_labels)
                    labels.extend(test_labels)

                    print("\nTrain Model")
                    print("==============\n")
                    print("Features: {} Labels: {}".format(features.shape, len(labels))) 
                    scaler = StandardScaler()
                    #grid = GridSearchCV(CLF, param_grid, refit=True, verbose=0)
                    #grid.fit(features, labels)
                    #print("Best Parameters")
                    #print(grid.best_params_)
                    
                    #Stack Models inside round
                    prediction1 = cross_val_predict(CLF1, scaler.fit_transform(features), labels, cv=15)
                    prediction2 = cross_val_predict(CLF2, scaler.fit_transform(features), labels, cv=15)
                    prediction3 = cross_val_predict(CLF3, scaler.fit_transform(features), labels, cv=15)
                    prediction4 = cross_val_predict(CLF4, scaler.fit_transform(features), labels, cv=15)
                    new_features = features.copy()
                    new_features["SVM"] = prediction1
                    new_features["linear"] = prediction2
                    new_features["KN"] = prediction3
                    new_features["GP"] = prediction4

                    
                    #Stack Models outside round
                    prediction_stacked = cross_val_predict(CLF5, scaler.fit_transform(new_features), labels, cv=15)
                    stack_features = features.copy() if not stack_models else stack_features
                    stack_features = pds.concat([stack_features, pds.DataFrame({"round_{}_svm".format(i):prediction_stacked})], axis=1) 
                    #stack_features = pds.concat([stack_features, pds.DataFrame({"round_{}_svm".format(i):prediction1, 
                    #     "round_{}_linear".format(i):prediction2, "round_{}_KN".format(i):prediction3,
                    #     "round_{}_GP".format(i): prediction4 } )], axis=1)
                    stack_models = True

                    #Final prediction:
                    from sklearn.impute import SimpleImputer
                    imputer = SimpleImputer(strategy="constant")
                    preprocess_data = scaler.fit_transform(imputer.fit_transform(stack_features))
                    prediction = cross_val_predict(CLF6, preprocess_data, labels, cv=15)
                    
                    #Score
                    value = median_absolute_error(np.ravel(labels), np.ravel(prediction))
                    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.ravel(labels), prediction)
                    print("Score:", value)
                    print("r2", r_value)
                    scores.append(value)
                    rsquares.append(r_value)
                    steps.append(j)

                    #Plot
                    fig, ax = plt.subplots()
                    ax.plot(steps, rsquares)
                    ax.plot(steps, scores)
                    fig.savefig("rsquare_scores.png")

                    fig, ax = plt.subplots()
                    ax.scatter(labels, prediction)
                    fig.savefig("result.png")

                    #Get results
                    results = np.absolute(np.ravel(prediction) - np.ravel(np.array(labels)))
                    
                    

                    #Reach convergenge
                    model1 = CLF1.fit(scaler.fit_transform(features), labels)
                    model2 = CLF2.fit(scaler.fit_transform(features), labels)
                    model3 = CLF3.fit(scaler.fit_transform(features), labels)
                    model4 = CLF4.fit(scaler.fit_transform(features), labels)
                    model_stacked = CLF5.fit(scaler.fit_transform(new_features), labels)
                    models.append([model1, model2, model3, model4, model_stacked])
                    previous_score = value
                    if r_value > 0.60:
                         preprocess_data = scaler.fit_transform(imputer.fit_transform(stack_features))
                         model_to_docking = CLF6.fit(scaler.fit_transform(preprocess_data), labels)
                         use_docking_model = True
                j+=1
        
    
def retrieve_features(df, Complex, molecules_3D, fragments, core, resname, positions_to_grow, bonded_atom, core_smile_pattern, labels=None, save=True):

   # Preprocess fragment
   features_test = pds.DataFrame(columns=COLUMNS) 
   true_labels = []
   print("nmolecules: {}\n  nfragments: {},\n nlabels: {}\n".format(len(molecules_3D), len(fragments), len(labels)))
   for i, (mol, fragment, bonded, label) in tqdm(enumerate(zip(molecules_3D, fragments, bonded_atom, labels)), total=len(fragments)):
       #Adapt ring to the right index
       next_ring = str(max([int(c) for c in core_smile_pattern if c.isdigit()]) + 1)
       for i in range(1, int(next_ring)):
           new_fragment = fragment.replace(str(i), next_ring)
           next_ring = str(int(next_ring) +  1)
           if new_fragment:
               fragment = new_fragment
           
       #mol = Chem.AddHs(mol)
       rec_info = []
       frag_info = [ 1 if i==bonded else 0 for i in range(positions_to_grow) ]
       atoms1 = pd.parsePDB(Complex)
       fragment = Chem.MolFromSmiles(fragment.strip("!"))

       # Count Heteroatoms
       n_carbon, n_nitrogen, n_oxigen = 0, 0, 0
       for atom in fragment.GetAtoms():
           if atom.GetAtomicNum() == 6: n_carbon += 1
           if atom.GetAtomicNum() == 7: n_nitrogen += 1
           if atom.GetAtomicNum() == 8: n_oxigen += 1
       frag_info.extend([n_carbon, n_nitrogen, n_oxigen])       

       # Fragment Properties
       fragment_features = collections.OrderedDict({ "Donor": 0, "Acceptor":0, "Aromatic": 0,  "Hydrophobe":0})
       for feature in factory.GetFeaturesForMol(fragment):
           feature_family  = feature.GetFamily()

       #Embed fragment properties
           if feature_family in fragment_features:
               fragment_features[feature_family] += 1
       # Calculate fragmetn volume
       #fragment_h = Chem.AddHs(fragment)
       #AllChem.EmbedMolecule(fragment_h)
       #fragment_volume = Chem.AllChem.ComputeMolVolume(fragment_h)
       #Embed information on a vector
       #frag_info = fragment_volume + fragment__features.values() 
       frag_info.extend(fragment_features.values())
       frag_info.append(len(fragment.GetAtoms()))
       frag_MACCS_fp = MACCSkeys.GenMACCSKeys(mol).ToBitString()
       for i in range(1, 167):
           frag_info.append(int(frag_MACCS_fp[i]))
       frag_info.append(Chem.rdMolDescriptors.CalcNumRotatableBonds(fragment))
       frag_info.append(Chem.rdMolDescriptors.CalcNumRings(fragment))
       frag_info.append(Chem.rdMolDescriptors.CalcNumHBD(fragment))
       frag_info.append(Chem.rdMolDescriptors.CalcNumHBA(fragment))
       frag_info.append(Chem.rdMolDescriptors.CalcNumHeteroatoms(fragment))
       frag_info.append(Chem.rdMolDescriptors.CalcNumAmideBonds(fragment))
       frag_info.append(Chem.rdMolDescriptors.CalcFractionCSP3(fragment))
       frag_info.append(Chem.rdMolDescriptors.CalcExactMolWt(fragment))
     
       # Align molecules
       molec = rdMolAlign.GetCrippenO3A(mol, core)
       molec.Align() 

       #Centroid of fragment
       fragment_atom_indxs = mol.GetSubstructMatch(fragment)
       molec_atoms = mol.GetAtoms()
       coords = [mol.GetConformer().GetAtomPosition(i) for i, atom in enumerate(molec_atoms) if i in fragment_atom_indxs ]
       try:
           fragment_centeroid = centeroidnp(np.array(coords))
       except IndexError:
           print("Molecule {} skipped".format(mol.GetProp("_Name")))
           continue
       frag_info = np.hstack([frag_info, fragment_centeroid])
       #Extract receptor info
       closer_residues = pd.parsePDB(Complex).select("protein within 10 of center", center=fragment_centeroid)
       distances = [ np.linalg.norm(fragment_centeroid-coords) for coords in closer_residues.getCoords() ]
       min_dist, average_dist, n_close_residues = np.min(distances), np.mean(distances), len(distances)
       residues = [0] * len(VOC_AMINOACIDS)
       for res in closer_residues.getResnames():
           idx = VOC_AMINOACIDS.index(res)
           residues[idx] += 1
       rec_info = np.hstack([np.array([min_dist, average_dist, n_close_residues]), np.array(residues)])
       fragment_receptor_info = np.hstack([frag_info, rec_info])
       
       # Save to general df
       df = df.append(pds.Series(fragment_receptor_info, index=COLUMNS), ignore_index=True)
       features_test = features_test.append(pds.Series(fragment_receptor_info, index=COLUMNS), ignore_index=True)
       true_labels.append(label)
   if save:
       df_save = features_test.copy()
       df_save["labels"] = true_labels
       df_save.to_csv("descriptors.csv")
   return df, features_test, true_labels

def centeroidnp(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return np.array([sum_x/length, sum_y/length, sum_z/length])
        

def dock(Complex, ligands):
    docking_obj = gl.Glide_Docker([Complex,], [ligands,]) 
    docking_obj.dock()
      

def combine_smiles(core_smiles_pattern, fragments, positions_to_grow, n_atoms):
    molecules = []
    bonded_atom = []
    molecules_fragment = []
    next_ring = str(max([int(c) for c in core_smiles_pattern if c.isdigit()]) + 1)
    for i in range(positions_to_grow):
        for fragment in fragments:
            if len(fragment.strip("!")) < n_atoms:
                continue
            for i in range(int(next_ring)):
                fragment = fragment.replace("1", next_ring)
            new_smiles = core_smiles_pattern.replace("*", "(" + fragment.strip("!") + ")", 1).replace("*", "")
            new_molecule = Chem.MolFromSmiles(new_smiles)
            if new_molecule and new_molecule not in molecules:
                molecules.append(new_molecule)
                bonded_atom.append(i)
                molecules_fragment.append(fragment)
        core_smiles_pattern = core_smiles_pattern.replace("*", "", 1)
    return molecules, molecules_fragment, bonded_atom
    
def identify_spots_to_grow(smiles):
    final_smiles = smiles[:]
    inserted = 0
    for i in range(len(smiles)):
        try:
            if smiles[i+1].isdigit():
                continue
        except IndexError:
            pass
        smiles_to_try = smiles[:]
        smiles_to_try = smiles[0:i+1] + "(C)" + smiles[i+1:]
        if Chem.MolFromSmiles(smiles_to_try):
            final_smiles = final_smiles[0:i+inserted+1] + "*" + final_smiles[i+inserted+1:] 
            inserted += 1
    print("Positions apt for growing: {}".format(inserted))
    print("Core smiles: {}".format(final_smiles))
    return final_smiles, inserted

def grow_molecule(smiles, vocabulary,  model, max_len=25, vocab_size=13):
    numerical_dataset = [vocabulary[c] for c in smiles]
    numerical_dataset = pad_sequences([numerical_dataset], maxlen=max_len)
    X_pred = np.array(numerical_dataset)
    X_pred_trans = to_categorical(X_pred, num_classes=vocab_size)
    valueToFind = model.predict(X_pred_trans.reshape(1, max_len, vocab_size))
    carac_index = np.argsort(valueToFind)[0][-4:]
    new_smiles = []
    for carac in carac_index:
      for item  in vocabulary.items():
          if item[1] == carac:
              new_smile = smiles + item[0]
              new_smiles.append(new_smile)
    return new_smiles

def retrieve_labels(molecules, dock_file):
    with open(dock_file, "r") as f:
        #Skip first line and report last column
        lines = f.readlines()[1:]

        #Get molecule indexes
        mols = []
        labels = []
        for line in lines:
            name = int(line.split(",")[-2].strip("molecule"))
            label = float(line.split(",")[-1].strip("\n"))
            if name not in mols:
                mols.append(name)
                labels.append(label)
        indexes_mols = np.array(mols)
        labels = np.array(labels)

        #Discard non dock molecules to learn
        final_labels = []
        final_molecules = []
        for mol in molecules:
            name = int(mol.GetProp("_Name").strip("molecule"))
            if name in indexes_mols:
                i, = np.where(indexes_mols == name)
                final_labels.append(labels[i])
                final_molecules.append(mol)

        return final_molecules, final_labels
  
def create_model(hidden_size = 800, dropout=0.2, batch_norm=1, lr=0.001, init_mode="uniform", activation="softplus", max_len=25, vocab_size=13):
   
  model = Sequential()
  model.add(Bidirectional(LSTM(hidden_size, kernel_initializer=init_mode), input_shape=(max_len, vocab_size)))
  model.add(Dropout(dropout))
  if batch_norm == 1:
    model.add(BatchNormalization())
  model.add(Dense(vocab_size, activation=activation))

  opt = RMSprop(lr=lr, clipnorm=1)
  model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])
  model.load_weights(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'model/weights/final_subset_cpu.hdf5'))
  print(model.summary())
  return model 
  
def parse_args(parser):
    parser.add_argument('--pdb',  type=str, help='Complex pdb of receptor + core to grow on')
    parser.add_argument('--resname', type=str, help='resname of the core')
    parser.add_argument('--only_grow', action="store_true", help='Perform only growin')
    parser.add_argument('--only_rank', action="store_true", help='Perform only ranking')
    parser.add_argument('--grow_iterations', type=int, help='Growing iteration (Same as maximum number of grown atoms)', default=10)

def main(pdb, resname, only_grow=True, only_rank=True, iterations=10):
    pdb = os.path.abspath(pdb)
    ligand_pdb, receptor_pdb = split_complex(pdb, resname)
    ligand_sdf = pdb_to_sdf(ligand_pdb)
    if only_grow:
        grow_protocol(pdb, ligand_sdf, VOC, resname, grow=True, rank=False, iterations=iterations)
    elif only_rank:
        grow_protocol(pdb, ligand_sdf, VOC, resname, grow=False, rank=True, iterations=iterations)
    else:
        grow_protocol(pdb, ligand_sdf, VOC, resname, grow=True, rank=True, iterations=iterations)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='GrowingAI')
    parse_args(parser)
    args = parser.parse_args()
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='GrowingAI')
    parse_args(parser)
    args = parser.parse_args()
    main(args.pdb, args.resname, args.only_grow, args.only_rank, args.grow_iterations)
