"""Main Script to run the MolClus pipeline
"""
import os
import argparse
from molclust.fingerprints import FingerprintCalculator
import pandas as pd
from molclust.clustering import Cluster

def main(dataframe:pd.DataFrame, 
        cluster:str,
         smilescol:str="SMILES", 
         n_clusters:int=None, 
         threshold:float=None,
         ) -> pd.DataFrame:
    """Main CLI

    Args:
        dataframe(pd.Dataframe): Dataframe containg  at least the SMILES strings 

    Returns:
        dataframe (pd.Dataframe): Dataframe with the labels column added
    """

    if cluster=='Butina':
        labels = Cluster(fpSize=1024, radius=2).Butina(smiles=dataframe[smilescol], threshold=threshold)
        cluster_map = {}
        for cluster_id, cluster in enumerate(labels, start=1):
            for row_idx in cluster:
                cluster_map[row_idx] = cluster_id
        dataframe['labels'] = dataframe.index.map(cluster_map)
    else:
        labels = Cluster(fpSize=1024, radius=2).Kmeans(dataframe[smilescol], n_clusters=n_clusters, normalize=True)
        dataframe['labels'] = labels

    print('Clutering done!')
    return dataframe 


if __name__=="__main__":
    p= argparse.ArgumentParser(description="Arguments for the clustering pipeline")
    p.add_argument('--dataframe', required=True, type=str, help='Path to the .csv file containing the molecules (strings) data')
    p.add_argument('--smilescol', type=str, default='SMILES', help='Name for the column containg the smiles')
    p.add_argument('--threshold', type=float,  help='Cuttof similarity distance to be used in Butina clustering. Can range from 0.1 to 1.0 :')
    p.add_argument('--cluster', type=str, required=True, choices={'Butina', 'KMeans'}, help='Clustering algorithm to be used')
    p.add_argument('--num-clusters', type=int, help='Number of clusters to form using KMeans algorithm')
    p.add_argument('--output', type=str, default='.', help='Output path for the new pandas dataframe')
    args = p.parse_args()

    
    assert args.dataframe.split('.')[-1] == 'csv', "Only .csv files are supported. Introduce the full path to a correct .csv file (e.g. `my_file.csv`)"

    try:
        datafame_name = args.dataframe.split('/')[-1].split('.')[0]
    except:
        dataframe_name = 'output_dataframe'
    try:
        dataframe = pd.read_csv(args.dataframe)
    except Exception as e:
        raise ValueError("An exception ocurred when reading the dataframe: {e}")
    print(f"Dataframe with {len(dataframe)} molecules loaded successfully")

    if args.cluster == 'Butina':
        assert args.threshold is not None, "When Butina algorithm is selected a threshold must be passed. Do so using `--threshold`"
        assert args.threshold > 0.1 and args.threshold < 1.0, f"The threshold for Butina algorithm must be passed and be between 0.1 to 1.0. Got {args.threshold}"
    else:
        assert args.num_clusters > 0, "The number of clusters must be a poositive number. Use `python main.py --dataframe my_df.csv --cluster KMeans --num-clusters 100`"
        assert args.num_clusters < len(dataframe), "There cannot be more clusters than molecules in the dataframe"



    main(dataframe,
         args.cluster,
         args.smilescol,  
         args.num_clusters, 
         args.threshold
         ).to_csv(os.path.join(args.output, f'{datafame_name}-clustered.csv'), index=False)