# Main class for clustering

import numpy as np 
from rdkit import Chem 
from tqdm import tqdm
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from typing import List, Optional
from molclust.fingerprints import FingerprintCalculator
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

class Cluster():

    def __init__(self, 
                 fpSize:int=1024, 
                 radius:int=2):
        self.fp_calc = FingerprintCalculator() 
        self.fpSize = fpSize
        self.radius = radius


    def Butina(self, 
               smiles:List[str], 
               threshold:float
               )-> List[int]:
        """Butina algorithm for clustering of molecules. Returns labels 

        Args:
            smiles (List[str]): List of the SMILES to be used in the clustering
            threshold (int): Cuttof similarity distance to be used in the clustering. Can range 
            from 0.1 to 1.0 
        Returns:
            labels (List[int]): Cluster label for every smiles passed 
        """

        assert threshold > 0.1 and threshold < 1.0, "The threshold should be between 0.1 and 1.0"

        matrix = []

        fp_arr = self.fp_calc.FingerprintFromSmiles(smiles, 
                                                    fp='morgan', 
                                                    fpSize = self.fpSize,
                                                    radius = self.radius, 
                                                    to_numpy= False)

        for i in tqdm(range(len(smiles)), desc='Clustering molecules', total=len(smiles)):
            for j in range(i):
                simil = DataStructs.FingerprintSimilarity(fp_arr[i], fp_arr[j])
                matrix.append(1 - simil)

        clusters= Butina.ClusterData(data=matrix, nPts=len(smiles), distThresh=(1 - threshold), isDistData=True)
        clusters= sorted(clusters, key=len, reverse=True)

        print(f"{len(set(clusters))} clusters formed")

        return clusters 



    def Kmeans(self, 
               smiles:List[str], 
               n_clusters: int,
               normalize:bool=True
               )-> List[int]:
        """KMeans algorithm for clustering of molecules. Returns labels 

        Args:
            smiles (List[str]): List of the SMILES to be used in the clustering
        Returns:
            labels (List[int]): Cluster label for every smiles passed 
        """

        fp_arr = self.fp_calc.FingerprintFromSmiles(smiles, 
                                                    fp='mqn', 
                                                    to_numpy= True)

        if normalize:
            fp_arr = StandardScaler().fit_transform(fp_arr)

        kmeans= KMeans(n_clusters=n_clusters)
        labels = kmeans.fit_predict(fp_arr) 

        print(f"{len(set(labels))} clusters formed")

        return labels 
        


