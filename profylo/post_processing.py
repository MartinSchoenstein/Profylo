import profylo.pre_processing as pp
import numpy as np
import networkx as nx
from Bio import Phylo
import seaborn as sns
import matplotlib.patches as mpatches
from sklearn.cluster import AgglomerativeClustering
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
import os
import urllib.request
from goatools.obo_parser import GODag
from goatools.associations import read_gaf
from goatools.test_data.genes_NCBI_9606_ProteinCoding import GENEID2NT
from goatools.goea.go_enrichment_ns import GOEnrichmentStudy
import markov_clustering as mc
from networkx.algorithms.community import fast_label_propagation_communities
from ete3 import Tree, TreeNode
import warnings
warnings.filterwarnings('ignore')


def _input_modules(x):
    if isinstance(x, list):
        return x
    elif isinstance(x, str):
        with open(x, "r") as f:
            liste = []
            for l in f:
                liste.append(l.strip().split(","))
        return liste
    else:
        ("Only accepted types for x are: a path to a txt file or a list type variable")

def _input_module(x):
    if isinstance(x, list):
        return x
    elif isinstance(x, str):
        with open(x, "r") as f:
            for l in f:
                liste = (f.strip().split(","))
        return liste
    else:
        ("Only accepted types for x are: a path to a txt file or a list type variable")

def hierarchical_dendrogram(
    x,                       
    method,                  
    path = None              
):
    """Function to obtain the dendrogram linked to hierarchical clustering from distances between profiles

    Args:
        x (str, pd.DataFrame): Distance matrix
        method (str): Method to use for hierarchical clustering : "ward", "average", "weighted", "complete"
        path (str, optional): Path to download the dendrogram plot. Defaults to None.
    """
    dfx = pp._input(x, test_binary = False)
    dfx = dfx.replace([np.inf, -np.inf], np.nan).fillna(0)
    dfx[dfx < 0] = 0
    dfx = 1 - dfx
    np.fill_diagonal(dfx.values, 0)
    dfx = squareform(dfx)
    L = linkage(dfx, method)
    plt.figure(figsize=(25, 10))
    dendrogram(L)
    plt.show
    if path is not None:
        plt.savefig(path) 


def hierarchical_clustering(
    x,                       
    method,                  
    criterion,               
    threshold,               
    path = None              
):
    """Function for performing hierarchical clustering on distances between profiles

    Args:
        x (str, pd.DataFrame): Distance matrix
        method (str): Method to use for hierarchical clustering : "ward", "average", "weighted", "complete"
        criterion (str): The way to cut the clusters, by "distance" between clusters or by the "number" of clusters 
        threshold (int): Threshold to use for cutting clusters
        path (str, optional): Path to download the clusters. Defaults to None.

    Returns:
        list: Returns a list containing the list of genes in each clusters
    """
    dfx =  pp._input(x, test_binary = False)
    dfx = dfx.apply(pd.to_numeric, errors="coerce").fillna(0)
    dfx = squareform(dfx)
    L = linkage(x, method)
    if criterion == "distance":
            clusters = fcluster(L, threshold, "distance")
    if criterion == "number":
            clusters = fcluster(L, threshold, "maxclust")
    clusters_number = len(np.unique(clusters))
    cluster_recap = [[] for _ in range(clusters_number)]
    names = x.columns.tolist()
    for i, cluster_id in enumerate(clusters):
        cluster_recap[cluster_id - 1].append(names[i])
    if path is not None:
        with open(path, "w") as f:
            for c in cluster_recap:
                f.write(",".join(c) + "\n")
    return cluster_recap


def graph_modules(
    x,                       
    distance,                
    threshold,               
    path = None              
):
    """Function to represent distances between profiles in a graph and extract connected components as modules

    Args:
        x (str, pd.DataFrame): Distance matrix
        distance (str): Metric used to obtain distances
        threshold (int): Trust threshold for creating edges
        path (str, optional): Path to download the modules. Defaults to None.

    Returns:
        list: Returns a list containing the list of genes in each modules
    """
    dfx = pp._input(x, test_binary = False)
    if distance in ["cotransition", "pearson", "mi", "Cotransition", "Pearson", "MI"] :
        dfx = dfx.replace([np.inf, -np.inf], np.nan).fillna(0)
        x_mask = dfx.mask(dfx < threshold, 0)
        np.fill_diagonal(x_mask.values, 0)
    if distance in ["svd_phy", "jaccard", "hamming", "SVD_phy", "Jaccard", "Hamming"]:
        dfx = dfx.replace([np.inf, -np.inf], np.nan).fillna(1)
        x_mask = dfx.mask(dfx > threshold, 1)
        np.fill_diagonal(x_mask.values, 1)
        x_mask = 1 - x_mask
    if distance in ["pcs", "PCS"]:
        dfx = dfx.replace([np.inf, -np.inf], np.nan).fillna(0)
        x_mask = dfx.mask(dfx < threshold, 0)
        max = x_mask.max()
        x_mask = x_mask/max
        np.fill_diagonal(x_mask.values, 0)
    G = nx.from_pandas_adjacency(x_mask)
    G.remove_edges_from([(u, v) for u, v, d in G.edges(data=True) if d['weight'] == 0])
    modules_recap = []
    for c in sorted(nx.connected_components(G), key=len, reverse=True):
        modules_recap.append(list(c))
    if path is not None:
        with open(path, "w") as f:
            for m in modules_recap:
                f.write(",".join(m) + "\n")
    return modules_recap


def markov_clustering(
    x,                       
    distance,                
    path = None              
):
    """Function to represent distances between profiles in a graph and performe markov clustering

    Args:
        x (str, pd.DataFrame): Distance matrix
        distance (str): Metric used to obtain distances
        path (str, optional): Path to download the modules. Defaults to None.

    Returns:
        list: Returns a list containing the list of genes in each modules
    """
    dfx = pp._input(x, test_binary = False)
    if distance in ["cotransition", "pearson", "mi", "Cotransition", "Pearson", "MI"] :
        dfx = dfx.replace([np.inf, -np.inf], np.nan).fillna(0)
        dfx[dfx < 0] = 0
        np.fill_diagonal(dfx.values, 0)
    if distance in ["svd_phy", "jaccard", "hamming", "SVD_phy", "Jaccard", "Hamming"]:
        dfx = dfx.replace([np.inf, -np.inf], np.nan).fillna(1)
        np.fill_diagonal(dfx.values, 1)
        dfx = 1 - dfx
    if distance in ["pcs", "PCS"]:
        dfx = dfx.replace([np.inf, -np.inf], np.nan).fillna(0)
        dfx[dfx < 0] = 0
        max = dfx.max()
        dfx = dfx/max
        np.fill_diagonal(dfx.values, 0)
    G = nx.from_pandas_adjacency(dfx)
    G.remove_edges_from([(u, v) for u, v, d in G.edges(data=True) if d['weight'] == 0])
    matrix = nx.to_numpy_array(G)
    result = mc.run_mcl(matrix)
    clusters = mc.get_clusters(result)
    genes = dfx.index.tolist()
    clusters_recap = []
    for c in clusters:
        liste = []
        for g in c:
            liste.append(genes[g])
        clusters_recap.append(liste)
    if path is not None:
        with open(path, "w") as f:
            for m in clusters_recap:
                f.write(",".join(m) + "\n")
    return clusters_recap


def label_propagation(
    x,                       
    distance,                
    threshold,               
    seed = None,             
    path = None              
):
    """Function to represent distances between profiles in a graph and performe label propagation clustering

    Args:
        x (str, pd.DataFrame): Distance matrix
        distance (str): Metric used to obtain distances
        threshold (int): Trust treshold to use for creating edges
        seed (int, optional): Seeding for the random part of the label propagation algorithm. Defaults to None.
        path (str, optional): Path to download the modules. Defaults to None.

    Returns:
        list: Returns a list containing the list of genes in each modules
    """
    dfx = pp._input(x, test_binary = False)
    if distance in ["cotransition", "pearson", "mi", "Cotransition", "Pearson", "MI"] :
        dfx = dfx.replace([np.inf, -np.inf], np.nan).fillna(0)
        x_mask = dfx.mask(dfx < threshold, 0)
        np.fill_diagonal(x_mask.values, 0)
    if distance in ["svd_phy", "jaccard", "hamming", "SVD_phy", "Jaccard", "Hamming"]:
        dfx = dfx.replace([np.inf, -np.inf], np.nan).fillna(1)
        x_mask = dfx.mask(dfx > threshold, 1)
        np.fill_diagonal(x_mask.values, 1)
        x_mask = 1 - x_mask
    if distance in ["pcs", "PCS"]:
        dfx = dfx.replace([np.inf, -np.inf], np.nan).fillna(0)
        x_mask = dfx.mask(dfx < threshold, 0)
        max = x_mask.max()
        x_mask = x_mask/max
        np.fill_diagonal(x_mask.values, 0)
    G = nx.from_pandas_adjacency(x_mask)
    G.remove_edges_from([(u, v) for u, v, d in G.edges(data=True) if d['weight'] == 0])
    F = fast_label_propagation_communities(G, seed = seed, weight = None)
    label_recap = []
    for c in sorted(F, key=len, reverse=True):
        label_recap.append(list(c))
    if path is not None:
        with open(path, "w") as f:
            for m in label_recap:
                f.write(",".join(m) + "\n")
    return label_recap


def go_enrichment(
    x,                       
    gaf,                 
    path = None,             
    complete_results = False 
):
    """Function to process GO enrichment test on a list of genes

    Args:
        x (str, list): List of list of genes (UNIPROT ID) or txt file with a line for each cluster
        gene2go (str): Link to a gene2go file like "http://current.geneontology.org/annotations/goa_human.gaf.gz"
        path (str, optional): Path to a dir to download full enrichment results. Defaults to None.
        complete_results (bool, optional): If =True, download full set of results files comming from goatools. Defaults to False.

    Returns:
        pd.DataFrame: Return resume results for each module
    """
    go_obo_url = "http://current.geneontology.org/ontology/go-basic.obo"
    if not os.path.exists("go.obo"):
        urllib.request.urlretrieve(go_obo_url, "go.obo")
    godag = GODag("go.obo")
    gene2go = read_gaf(gaf, godag=godag, namespace=None)
    genes_fond = set(gene2go.keys())
    x = _input_modules(x)
    df = pd.DataFrame(columns = ["length","GO/P-value(bh)_1", "GO/P-value(bh)_2", "GO/P-value(bh)_3", "GO/P-value(bh)_4", "GO/P-value(bh)_5"])
    for i, module in enumerate(x):
        goeaobj = GOEnrichmentStudy(genes_fond, gene2go, godag, propagate_counts=True, alpha=0.05, methods=['fdr_bh'])
        results = goeaobj.run_study(module)
        results = sorted(results, key = lambda x: x.p_fdr_bh)[:5]
        data = {}
        for j, res in enumerate(results):
            data["length"] = len(module)
            data["GO/P-value(bh)_" + str(j+1)] = [res.GO, res.name, res.p_fdr_bh]
        df.loc[len(df)] = data
        if complete_results is True:
            path2 = path + "/" + str(i) + ".tsv"
            goeaobj.wr_tsv(path2, results)
    if path is not None and complete_results is True:
        path1 = path + "/resume_results.csv"
        df.to_csv(path1, index = False)
    if path is not None and complete_results is False:
        df.to_csv(path, index = False)
    return df


def profils_heatmap(
    x,                       
    selection = None,        
    tree = None,             
    clades = None,           
    ordered = False,         
    path = None              
):
    """Function to represent profiles (presence/absence) on a heatmap

    Args:
        x (str, pd.DataFrame): Profile matrix
        selection (list, optional): Genes selection to represent. Defaults to None.
        tree (str, optional): Newick tree to order profils, not necessary if ordered = True. Defaults to None.
        clades (list, optional): Clades to highlight on heatmaps. Defaults to None.
        ordered (bool, optional): True if profils are already ordered else need a Newick tree. Defaults to False.
        path (str, optional): Path to download the heatmap. Defaults to None.
    """
    if ordered == False and tree == None:
        raise ValueError("Need a tree to order profiles")
    if clades is not None and tree == None:
        raise ValueError("Need tree to obtain phylogenetic informations")
    profils = pp._input(x, test_binary = False)
    if ordered == False :
        profils = pp.order_by_tree(profils, tree)
    if selection is None:
        profils_selection = profils
    else:
        profils_selection = profils.loc[selection]
    profils_selection = profils_selection.fillna(0)
    profils_selection = profils_selection.apply(pd.to_numeric)
    tree = Phylo.read(tree, "newick")
    if clades is not None:
        new_col = []
        for c in profils_selection.columns:
            add = False
            phylogenie = [clade.name for clade in tree.get_path(c) if clade.name]
            for clade in phylogenie:
                if str(clade) in clades:
                    if add == False:
                        new_col.append(clade)
                        add = True
            if add == False:
                new_col.append("other")
        df_clades = pd.DataFrame([new_col], columns=profils_selection.columns, index = ["clade"])
        cate = df_clades.loc["clade"]
        unique_categories = cate.unique()
        palette = dict(zip(unique_categories, sns.color_palette("Paired", len(unique_categories))))
        col_colors = cate.map(palette)
        handles = [mpatches.Patch(color=color, label=category) for category, color in palette.items()]
    plt.figure(figsize=(12, 6))
    if clades is not None:
        sns.clustermap(profils_selection, cmap="coolwarm", row_cluster=False, col_cluster=False, col_colors=col_colors, xticklabels=False, cbar_pos=None, dendrogram_ratio=(.01, .1))
        plt.legend(handles=handles, title="Clades", bbox_to_anchor=(0.5, 1.05), loc='lower right', ncol = 3)
    else:
        sns.clustermap(profils_selection, cmap="coolwarm", row_cluster=False, col_cluster=False, xticklabels=False, cbar_pos=None, dendrogram_ratio=(.01, .1))
    if path is not None:
        plt.savefig(path, bbox_inches='tight')
    plt.show()


def _state_on_nodes(x):
    if hasattr(x, "state"):
        return(x.state)
    count = 0
    for child in x.children:
        state = _state_on_nodes(child)
        if state == 1:
            count = count + 1
    if count > 0:
        x.add_feature("state",  1)
    else:
        x.add_feature("state", 0)
    return(x.state)


def tree_annotation(
    x,                  
    profils,            
    path_tree,          
    path = None         
):
    """Function to reconstruct the evolutionary history of genes on a tree from profiles

    Args:
        x (list): List of genes or txt file with all genes in one line
        profils (str, pd.DataFrame): profile matrix
        path_tree (str): Path to Newick tree to use
        path (str, optional): Path to download the tree. Defaults to None.

    Returns:
        Tree: Annotated tree with "state" on nodes representing mean presence/absence in the cluster
    """
    profils = pp._input(profils, test_binary = False)
    x = _input_module(x)
    liste_tree = []
    for gene in x:
        unique = False
        t = Tree(path_tree, format = 8)
        for leaf in TreeNode.iter_leaves(t):
            leaf.add_feature("state", profils.loc[gene, leaf.name])
        AC = TreeNode.get_tree_root(t)
        _state_on_nodes(AC)
        if sum(profils.loc[gene]) == 1:
            unique = True
            oldest_node = ""
            all_kids = []
        else:
            higher_kids = 0
            for node in t.traverse():
                if node.state == 1:
                    count = 0
                    for child in node.children:
                        if child.state == 1:
                            count = count + 1
                    if count >= 2 :
                        if len(node.get_leaves()) > higher_kids:
                            higher_kids = len(node.get_leaves())
                            oldest_node = node
            all_kids = list(oldest_node.iter_descendants())
        for node in t.traverse():
            if node != oldest_node:
                if node.state == 1:
                    if unique is True:
                        if not node.is_leaf():
                            node.state = 0
                    elif node not in all_kids:
                        node.state = 0
        liste_tree.append(t)
    t_mean = Tree(path_tree, format = 8)
    for node in t_mean.traverse(strategy="postorder"):
        mean = 0
        for tree in liste_tree:
            n = tree&node.name
            mean = mean + int(n.state)
        mean = mean / len(liste_tree)
        node.add_feature("state", mean)
    if path is not None:
        t_mean.write(outfile=path, format=8, features=["state"])
    return t_mean


def phylogenetic_statistics(
    x,                  
    profils = None,            
    path_tree = None,          
    path = None,        
    dl_tree = False 
):
    """Function to compute statistics from a cluster by a phylogenetic tree

    Args:
        x (list, str): List of list of genes or txt file with a line for each cluster, or list of annotated Trees
        profils (str, pd.DataFrame, optional): Profile matrix. Defaults to None (Mandatory if x is not a list of Tree).
        path_tree (str, optional): Path to Newick tree to use. Defaults to None (Mandatory if x is not a list of Tree).
        path (str, optional): Path to download output dataframe, must be a dir if dl_tree is True. Defaults to None.
        dl_tree (bool, optional): If dl_tree is true, return the mean tree used for each cluster. Defaults to False.

    Returns:
        pd.DataFrame: Length, Parsimony, Number of Leaf with a presence, Last common ancestor
    """
    profils = pp._input(profils, test_binary = False)
    x = _input_modules(x)
    df = pd.DataFrame(columns=["Cluster", "Length", "Parsimony", "Presence", "LCA"])
    for i, module in enumerate(x):
        if isinstance(module, list):
            t_mean = tree_annotation(module, profils = profils, path_tree = path_tree)
        else:
            t_mean = module
        parsimony = 0
        presence = 0
        for node in t_mean.traverse(strategy = "postorder"):
            if node.state >= 0.5:
                if node.is_leaf():
                    presence = presence + 1
                LCA = node.name
                for child in node.children:
                    if child.state < 0.5:
                        parsimony = parsimony + 1
        if dl_tree is True and path is not None:
            path_tree_dl = path + "/" + str(i + 1) + ".nhx"
            t_mean.write(outfile=path_tree_dl, format=8, features=["state"])
        df.loc[i] = [i + 1, 1 if isinstance(module, Tree) else len(module), parsimony, presence, LCA]
    if path is not None:
        if dl_tree is True:
            path_df = path + "/phylogenetic_statistics.csv"
        else:
            path_df = path
        df.to_csv(path_df, index = False)
    return df