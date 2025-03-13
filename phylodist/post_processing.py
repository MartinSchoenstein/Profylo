import PhyloDist.pre_processing as pp
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
import os
import urllib.request
from goatools.obo_parser import GODag
from goatools.associations import read_gaf
from goatools.test_data.genes_NCBI_9606_ProteinCoding import GENEID2NT
from goatools.goea.go_enrichment_ns import GOEnrichmentStudy
import warnings
warnings.filterwarnings('ignore')


def input_modules(x):
    if isinstance(x, list):
        return x
    elif isinstance(x, str):
        with open(x, "r", encoding="utf-8") as f:
            l = l.strip().strip("{}")
            liste = [elem.strip("'") for elem in l.split(", ")]
        return liste
    else:
        ("Only accepted types for x are: a path to a graph_modules() output file or a list type variable")



#Function to obtain a dendrogram from distances for hierarchical clustering 
def hierarchical_dendrogram(
    x,                       #Distances to use
    method,                  #Method to use for hierarchical clustering : "ward", "average", "weighted", "complete"
    path = None              #Path to download the dendrogram plot
):
    dfx = pp.input(x, test_binary = False)
    L = linkage(dfx, method)
    plt.figure(figsize=(25, 10))
    dendrogram(L)
    plt.show
    if path is not None:
        plt.savefig(path) 


#Function for performing hierarchical clustering on distances, returns a list containing the list of genes in each clusters
def hierarchical_clustering(
    x,                       #Distance to use
    method,                  #Method to use for hierarchical clustering : "ward", "average", "weighted", "complete"
    criterion,               #The way to cut the clusters, by distance or by the number of clusters
    threshold,               #Threshold to use for cutting clusters
    path = None              #Path to download the clusters
):
    dfx =  pp.input(x, test_binary = False)
    dfx = dfx.apply(pd.to_numeric, errors="coerce").fillna(0)
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


#Function to represent distance in a graph and extract connected components as modules, returns a list containing the list of genes in each modules
def graph_modules(
    x,                       #Distance to use
    distance,                #Metric used to obtain distances
    threshold,               #Trust treshold to use for creating nodes
    path = None              #Path to download the modules
):
    dfx = pp.input(x, test_binary = False)
    if distance in ["cotransition", "pearson", "mi", "Cotransition", "Pearson", "MI"] :
        x_mask = dfx.mask(dfx < threshold, 0)
        np.fill_diagonal(x_mask.values, 0)
    if distance in ["svd_phy", "jaccard", "hamming", "SVD_phy", "Jaccard", "Hamming"]:
        x_mask = dfx.mask(dfx > threshold, 1)
        np.fill_diagonal(x_mask.values, 1)
        x_mask = 1 - x_mask
    if distance in ["pcs", "PCS"]:
        x_mask = dfx.mask(dfx < threshold, 0)
        max = x_mask.max()
        x_mask = x_mask/max
        np.fill_diagonal(x_mask.values, 0)
    G = nx.from_pandas_adjacency(x_mask)
    modules_recap = []
    for c in sorted(nx.connected_components(G), key=len, reverse=True):
        modules_recap.append(list(c))
    if path is not None:
        with open(path, "w") as f:
            for m in modules_recap:
                f.write(",".join(m) + "\n")
    return modules_recap


#Function to process GO enrichment test on a list of genes, return list of GO therms with their corrected p-values ​​when < 0.05, can dl full results with path
def go_enrichment(
    x,                       #List of genes (UNIPROT ID)
    path = None              #Path to download full enrichment results
):
    go_obo_url = "http://purl.obolibrary.org/obo/go.obo"
    go_gaf_url = "http://current.geneontology.org/annotations/goa_human.gaf.gz"
    if not os.path.exists("go.obo"):
        urllib.request.urlretrieve(go_obo_url, "go.obo")
    if not os.path.exists("goa_human.gaf"):
        urllib.request.urlretrieve(go_gaf_url, "goa_human.gaf.gz")
        os.system("gunzip goa_human.gaf.gz")
    godag = GODag("go.obo")
    gene2go = read_gaf("goa_human.gaf", godag=godag, namespace=None)
    genes_fond = set(gene2go.keys())
    x = input_modules(x)
    goeaobj = GOEnrichmentStudy(genes_fond, gene2go, godag, propagate_counts=True, alpha=0.05, methods=['fdr_bh'])
    results = goeaobj.run_study(x)
    data = [(res.GO, res.p_fdr_bh) for res in results if res.p_fdr_bh < 0.05]
    df = pd.DataFrame(data, columns = ["GO", "p-value(bh)"])
    if path is not None:
        goeaobj.wr_tsv(path, results)
    return df


#Function to represent profiles on a heatmap
def profils_heatmap(
    x,                       #Profiles to use
    selection = None,        #Genes whose profiles will be represented, if None all profils are represented
    tree = None,             #Newick tree to order profils, not necessary if ordered = True
    clades = None,           #Clades to highlight on heatmaps
    ordered = False,         #If profils are already ordered --> True, else need a Newick tree
    path = None              #Path to download the heatmap
):
    if ordered == False and tree == None:
        raise ValueError("Need a tree to order profiles")
    if clades is not None and tree == None:
        raise ValueError("Need tree to obtain phylogenetic informations")
    profils = pp.input(x, test_binary = False)
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
    plt.show()
    if path is not None:
        plt.savefig(path)