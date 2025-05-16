import pandas as pd
import numpy as np
from Bio import Phylo
import warnings
warnings.filterwarnings('ignore')



def _is_binary(df):
    binary = df.map(lambda x: x in [-1, 0, 1]).all().all()
    return binary


def _input(x, test_binary = True):
    if isinstance(x, str) :
        dfx = pd.read_csv(x, sep=",", index_col=0)
    elif isinstance(x, pd.DataFrame):
        dfx = x
    else:
        raise ValueError("Only accepted types for x are: a path to a csv file or a pandas dataframe")
    if test_binary == True:
        binary = _is_binary(dfx)
        return dfx, binary
    else:
        return dfx


def to_binary(
    x,                   
    threshold=0.5
):
    """Function to binarize continuous profiles

    Args:
        x (str, pd.DataFrame): Profiles matrix 
        threshold (float, optional): Tresholds to apply, >tresholds --> 1, <tresholds --> 0. Defaults to 0.5.

    Returns:
        pd.DataFrame: Returns a binary profiles matrix
    """
    dfx= input(x, test_binary = False)
    return (dfx > threshold).astype(int)


def normalize(
    x                    
):
    """Function to normalize continuous profiles

    Args:
        x (str, pd.DataFrame): Profiles matrix

    Returns:
        pd.DataFrame: Returns normalized continuous profiles matrix
    """
    dfx = input(x, test_binary = False)
    for i in dfx.index:
        dfx.loc[i] = dfx.loc[i]/max(dfx.loc[i])
    return dfx


def transition_vector(
    x,                   
    path = None,         
    from_outside = True
):
    """Function to convert profiles matrix into transition vectors matrix : 00110 --> 0010-1. 

    Args:
        x (str, pd.DataFrame): Ordered profiles matrix
        path (str, optional): Path to use to download transition vectors. Defaults to None.

    Returns:
        pd.DataFrame: Returns transition vectors matrix
    """
    if from_outside is True:
        dfx, binary_x = input(x)
        if binary_x is False:
           dfx = to_binary(dfx)
           print("Profiles were not binary, to_binary() was applied with 0.5 for threshold")
    else:
        dfx = x
    tv = np.zeros((len(dfx.index),(len(dfx.columns))))
    for j in range(0, len(dfx)):
        vec = dfx.iloc[j]
        pos_ori = np.asarray(vec[:-1])
        pos_new = np.asarray(vec[1:])
        trans_vec = pos_new-pos_ori
        tv[j] = np.pad(trans_vec, (1,0))
    tv = pd.DataFrame(tv, index=dfx.index, columns=dfx.columns)
    if path is not None:
        tv.to_csv(path, index = True)
    return tv


def order_by_tree(
    x,                    
    tree,                 
    path =  None          
):
    """Function to order profiles matrix with a Newick tree

    Args:
        x (str, pd.DataFrame): Profiles matrix
        tree (str): Newick tree used for ordering profiles
        path (str, optional): Path to use to download profiles. Defaults to None.

    Returns:
        pd.DataFrame: Returns a tree ordered profiles matrix
    """
    if tree is None:
        raise TypeError("A newick file is missing")
    dfx = input(x, test_binary = False)
    phylo = Phylo.read(tree, "newick")
    leaf = phylo.get_terminals()
    ordered_df = pd.DataFrame()
    for l in leaf:
        ordered_df[str(l.name)] = dfx[str(l.name)]
    if path is not None:
        ordered_df.to_csv(path, index = True)
    return ordered_df