from rdkit.SimDivFilters import rdSimDivPickers  
from tqdm import tqdm 
from rdkit.Chem import AllChem
import numpy as np
import umap.umap_ as umap
from rdkit import RDLogger                                                                                                                                                               
RDLogger.DisableLog('rdApp.*')
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.SaltRemover import SaltRemover
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

def santize_molecules(paths):
    """Sanitize molecules 

    Args:
        paths (list[str]): List with the paths to .sdf files

    Returns:
        dataframe (pd.DataFrame): Dataframe with sanitize and named molecules
    """
    # prepare santizers 
    salt_remover = SaltRemover() 
    uncharger    = rdMolStandardize.Uncharger()
    clean_data = []
    for path in paths:
        supplier = Chem.SDMolSupplier(path, sanitize=True, removeHs=False)
        for mol in tqdm(supplier, desc='Cleaning'):
            if mol is None:
                continue
            try:
                # sanitize
                Chem.SanitizeMol(mol, catchErrors=True)
                # rm salt 
                mol = salt_remover.StripMol(mol)
                # Neutralize
                mol = uncharger.uncharge(mol, )
                # Sanitize again
                Chem.SanitizeMol(mol)
                # canonical smile
                smi = Chem.MolToSmiles(mol, isomericSmiles=True)
                # Get name if exists
                name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
                clean_data.append((smi, name))
            except Exception:
                # discard exceptions
                continue

    df_clean = pd.DataFrame(clean_data, columns=["SMILES", "Name"])
    return df_clean


def _mols_and_fps(df, smiles_col="SMILES", radius=2, nBits=2048) -> pd.DataFrame:
    """Convert a SMILES column to RDKit mols + Morgan/ECFP4 fingerprints.
            ) -> (rdkit.mol, fingerprints, idx)

    Args:
        df(pd.Dataframe): Dataframe containg  at least the SMILES strings 
        smiles_col(str): Name of the column containing the smiles strings
        radius(int): Radius for ECFP fingerprints. Defaults to 2
        nBits(int): Size of the ECFP fingerprint. Default to 2048

    Returns:
        dataframe (pd.Dataframe): Dataframe with the labels column added
    """ 
    fps = []
    for i, smi in tqdm(enumerate(df[smiles_col]), desc='Calculating fingerprints', total=len(df)):
        try:
            m = Chem.MolFromSmiles(smi)
        except Exception as e:
            print(f'An exception ocurred with {smi}: {e}')
        if m:
            fps.append(AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits))
    return fps


def picker(n_pick:int, df:pd.DataFrame, smiles_col:str="SMILES", radius:int=2, nBits:int=2048, plot:bool=True, seed:int=42):
    """Select N molecules from a dataset using MaxMinPicker selection. Optinally visualize through TMAP

    Args:
        n (int): Number of molecules to select from the dataset
        df (pd.DataFrame): _description_
        smiles_col (str, optional): _description_. Defaults to "SMILES".
        radius (int, optional): _description_. Defaults to 2.
        nBits (int, optional): _description_. Defaults to 2048.
        plot (bool, optional): _description_. Defaults to True.
        seed (int, optional): _description_. Defaults to 42.
    """
    df = pd.read_csv('cleaned_molecules.csv')
    assert n_pick < len(df), f"The number of molecules to select {n} cannot be larger than the dataframe {len(df)}"
    fps= _mols_and_fps(df=df, smiles_col=smiles_col, radius=radius, nBits=nBits) 
    picker = rdSimDivPickers.MaxMinPicker()               
    picks = picker.LazyBitVectorPick(fps,           # list[ExplicitBitVect] (fingerprints)
                                    len(fps),      # pool size
                                    n_pick,        # how many to pick
                                    [],            # no initial seed set
                                    seed)            # RNG seed


    if plot:
        #convert ot binary again 
        arr = np.asarray([np.frombuffer(fp.ToBitString().encode(), 'S1').astype(np.int8)
                        for fp in fps])

        reducer = umap.UMAP(n_neighbors=15,      
                            min_dist=0.1,
                            metric='jaccard',
                            random_state=seed)

        coords = reducer.fit_transform(arr)     

        mask = np.zeros(len(fps), dtype=bool)
        mask[picks] = True          # True for the chosen molecules

        plt.figure(figsize=(20, 20), dpi=110) # If figure is to large change (20,20) to smaller value e.g(6,5)

        plt.scatter(coords[~mask, 0], coords[~mask, 1],
                    s=10, c='black', alpha=0.9, linewidths=0)

        # highlight: MaxMin subset in red bc why not
        plt.scatter(coords[mask, 0], coords[mask, 1],
                    s=15, c='red', alpha=0.7, linewidths=0)

        plt.title(f'UMAP projection of {len(df)} compounds\nred = MaxMin subset')
        plt.axis('off')             # chemical space has no intrinsic axes
        plt.tight_layout()
        plt.show()

    return df.iloc[picks]

if __name__=="__main__":
    paths = ['Innopharm_chemical_library.sdf', 'Welab_chemical_library.sdf']
    df = santize_molecules(paths)
    df = picker(n_pick=5000, df=df, smiles_col="SMILES", radius=2, nBits=2048, plot=True, seed=49)
    df.to_csv('selection.csv', index=False)