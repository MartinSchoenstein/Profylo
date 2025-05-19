import profylo.distances as dst
import profylo.pre_processing as pp
import profylo.post_processing as post
import warnings
import argparse
warnings.filterwarnings('ignore')

import argparse


def distance_profiles(
    method,                 
    x,                      
    y=None,                 
    type = "matrix",        
    confidence=1.5,         
    penalty=0.6,            
    truncation= 0.5,        
    consecutive = True,     
    tree = None,            
    path = None             
):
    """Main function of the library, allows access to 7 metrics for comparing phylogenetic profiles

    Args:
        method (str): Name of the distance method to use : "Jaccard", "Hamming", "Pearson", "MI", "PCS", "Cotransition", "SVD_phy"
        x (str, pd.DataFrame): Data to use
        y (str, pd.DataFrame, optional): Second date source if you want to compare the x profils to y profils and not all x profils to all x profils. Defaults to None.
        type (str, optional): Type of data source, classic profils matrix ("matrix") or transition vectors matrix ("transition_vector"). Defaults to "matrix".
        confidence (float, optional): PCS parameter, influence the positive weight of a double match. Defaults to 1.5.
        penalty (float, optional): PCS parameter, influence the negative weight of a double missmatch. Defaults to 0.6.
        truncation (float, optional): SVD-phy parameter, influence the data reduction. Defaults to 0.5.
        consecutive (bool, optional): Cotransition score parameter, when consecuive=False: only counts one transition in 2 directly consecutive ones. Defaults to True.
        tree (str, optional): Newick tree to order profils for PCS or Cotransition score. Defaults to None.
        path (str, optional): Path to download distance dataframe. Defaults to None.

    Returns:
        pd.DataFrame: Returns a distance matrix in a dataframe

    """
    if method not in ["Jaccard", 'jaccard', "Hamming", 'hamming', "Pearson", 'pearson', "MI", 'mi', "Cotransition", 'cotransition', "PCS", 'pcs', "SVD_phy", 'svd_phy']:
        raise ValueError("Method not accepted")
    if method in ["Cotransition", 'cotransition', "PCS", 'pcs'] and tree == None:
        if type not in ["Transition_vector", "Transition vector", "transition_vector", "transition vector"]:
            raise ValueError("Tree needed for this method")
    if type not in ["matrix", "Matrix", "Transition_vector", "Transition vector", "transition_vector", "transition vector"]:
        raise ValueError("type must be matrix or transition_vector")
    if method in ["SVD_phy", 'svd_phy'] and y != None:
        raise ValueError("SVD_phy only accepts 1 matrix")
    dfx, binary_x = pp._input(x)
    if y is None:
        dfy = None
        print(len(dfx), " profiles loaded into ", method, " distance.")
    else:
        dfy = pp._input(y, test_binary = False)
        print(len(dfx) + len(dfy), " profiles loaded into ", method, " distance.")
    if method in ["Jaccard", 'jaccard', "Hamming", 'hamming', "Cotransition", 'cotransition', "PCS", 'pcs'] and binary_x == False:
        dfx = pp.to_binary(dfx)
        if dfy is not None:
            dfy = pp.to_binary(dfy)
        print("Profiles were not binary, to_binary() was applied with 0.5 for threshold")
    print("Running...")
    if method == "Jaccard" or method == "jaccard":
        result = dst.jaccard(dfx, dfy)
    if method == "Hamming" or method == "hamming":
        result = dst.hamming(dfx, dfy) 
    if method == "Pearson" or method == "pearson":
        result = dst.pearson(dfx, dfy)
    if method == "MI" or method == "mi":
        result = dst.mi(dfx, dfy)
    if method == "cotransition" or method == "Cotransition":
        if type == "matrix":
            ordered_dfx = pp.order_by_tree(dfx, tree)
            tvx = pp.transition_vector(ordered_dfx, from_outside = False)
            if dfy is None :
                ordered_dfy = None
                tvy = None
            else :
                ordered_dfy = pp.order_by_tree(dfy, tree)
                tvy = pp.transition_vector(ordered_dfy, from_outside = False)
        else:
            tvx = dfx
            tvy = None
        result, p_value = dst.cotransition(tvx, tvy, consecutive)
    if method == "pcs" or method == "PCS":
        if type == "matrix":
            ordered_dfx = pp.order_by_tree(dfx, tree)
            tvx = pp.transition_vector(ordered_dfx, from_outside = False)
            if dfy is None :
                ordered_dfy = None
                tvy = None
            else :
                ordered_dfy = pp.order_by_tree(dfy, tree)
                tvy = pp.transition_vector(ordered_dfy, from_outside = False)
        else:
            tvx = dfx
            tvy = None
        result = dst.pcs(tvx, tvy, confidence, penalty)
    if method == "svd_phy" or method == "SVD_phy":
        result = dst.SVD_phy(dfx, truncation)
    if path is not None:
        if method == "cotransition" or method == "Cotransition":
            p_value.to_csv(path, index = True)
        result.to_csv(path, index=True),
    print("Done.")
    return result

def make_modules(x, clustering, method = None, criterion = None, threshold = None, distance = None, seed = None, path = None):
    if clustering == "label_propagation":
        if distance is None:
            raise ValueError("Distance metric used is requested")
        if threshold is None :
            raise ValueError("Treshold to build edges is requested")
        post.label_propagation(x, distance, threshold, seed, path)
    if clustering == "markov_clustering":
        if distance is None :
            raise ValueError("Distance metric used is requested")
        post.markov_clustering(x, distance, path)
    if clustering == "connected_components" or clustering == "graph_modules":
        if distance is None :
            raise ValueError("Distance metric used is requested")
        if threshold is None :
            raise ValueError("Threshold to build edges is requested")
        post.graph_modules(x, distance, threshold, path)
    if clustering == "hierarchical_clustering":
        if distance is None:
            raise ValueError("Distance metric used is requested")
        if criterion is None :
            raise ValueError("Criterion to cut clusters is requested")
        if threshold is None :
            raise ValueError("Threshold to cut clusters is requested")
        post.hierarchical_clustering(x, distance, criterion, threshold, method, path)


def phylogenetic_statistics(x, profils = None, path_tree = None, path = None, dl_tree = False):
    post.phylogenetic_statistics(x, profils, path_tree, path, dl_tree)

def parse_args():
    parser = argparse.ArgumentParser(description='Profylo: Phylogenetic profile distance manipulaion')
    subparsers = parser.add_subparsers(dest='mode', required=True, help='Mode selection')
    distance_parser = subparsers.add_parser('distance_compute', help='Compute distance between profiles')
    distance_parser.add_argument('-m', '--method', type=str, required=True, choices= ["Jaccard", 'jaccard', "Hamming", 'hamming', "Pearson", 'pearson', "MI", 'mi', "Cotransition", 'cotransition', "PCS", 'pcs', "SVD_phy", 'svd_phy'], help='Distance method to use: Jaccard, Hamming, Pearson, MI, PCS, Cotransition, SVD_phy')
    distance_parser.add_argument('-x', '--matrix', type=str, required=True, help='Path to the first profile matrix')
    distance_parser.add_argument('-o', '--output', type=str,  required=True, help='Path to save the output distance matrix')

    distance_parser.add_argument('-x2', '--matrix2', type=str, help='Path to the second profile matrix (Optionnal for cross-profile comparisons)')
    distance_parser.add_argument('-c', '--confidence', type=float, default=1.5, help='Optionnal - PCS parameter: confidence')
    distance_parser.add_argument('-p', '--penalty', type=float, default=0.6, help='Optionnal - PCS parameter: penalty')
    distance_parser.add_argument('-tr', '--truncation', type=float, default=0.5, help='Optionnal - SVD_phy parameter: truncation')  
    distance_parser.add_argument('-co', '--consecutive',  action="store_true", help='Optionnal -Cotransition parameter: consecutive')
    distance_parser.add_argument('-tree', '--tree', type=str, help='Required for cotransition and PCS - Path to the newick tree file')

    modules_parser = subparsers.add_parser('make_modules', help='Make modules from a distance matrix')
    modules_parser.add_argument('-x', '--matrix', type=str, required=True, help='Path to the distance matrix')  
    modules_parser.add_argument('-cl', '--clustering', type=str, required=True, choices=["label_propagation", "markov_clustering", "connected_components", "graph_modules", "hierarchical_clustering"], help='Clustering method to use: label_propagation, markov_clustering, connected_components, graph_modules, hierarchical_clustering')
    modules_parser.add_argument('-o', '--output', type=str, required=True, help='Path to save the output modules')
    modules_parser.add_argument('-d', '--distance', type=str, required=True, help='Distance metric used for build the distance matrix.')
    modules_parser.add_argument('-m', '--method', type=str, default="ward", help='Only if method is hierarchical_clustering. Method to use for hierarchical clustering: single, complete, average, weighted, centroid, median, ward (default)')
    modules_parser.add_argument('-cr', '--criterion', type=str, help='Only if method is hierarchical_clustering. Criterion to use for hierarchical clustering: distance, maxclust')
    modules_parser.add_argument('-th', '--threshold', type=float, help='Threshold to use for graph construction - distance above the treshold or similarity below the treshold will not be included as edge.')
    modules_parser.add_argument('-s', '--seed', type=int, help='Label propagation only. Seed to use for label propagation.')
    modules_parser = subparsers.add_parser('phylogenetic_statistics', help='Compute phylogenetic statistics')
    modules_parser.add_argument('-m', '--modules', type=str, required=True, help='Text file of list of genes, with a line for each cluster')
    modules_parser.add_argument('-p', '--profil', type=str, required=True, help='Path to the profile matrix')
    modules_parser.add_argument('-tree', '--path_tree', required=True, type=str, help='Path to the newick tree file')
    modules_parser.add_argument('-o', '--output', required=True, type=str, help='Path to save the output phylogenetic statistics')
    modules_parser.add_argument('-dl', '--dl_tree',   action="store_true", help='Use if this function should output an annotated by module')

    args = parser.parse_args()
    return args

def profylo_cli():
    args = parse_args()

    if args.mode == 'distance_compute':
        distance_profiles(
            method=args.method,
            x=args.matrix,
            y=args.matrix2,
            type="matrix",
            confidence=args.confidence,
            penalty=args.penalty,
            truncation=args.truncation,
            consecutive=args.consecutive,
            tree=args.tree,
            path=args.output)
    elif args.mode == 'make_modules':
        make_modules(
            x=args.matrix,
            clustering=args.clustering,
            method=args.method,
            criterion=args.criterion,
            threshold=args.threshold,
            distance=args.distance,
            seed=args.seed,
            path=args.output)
    elif args.mode == 'phylogenetic_statistics':
        phylogenetic_statistics(
            x=args.modules,
            profils=args.profil,
            path_tree=args.path_tree,
            path=args.output,
            dl_tree=args.dl_tree)   
def _phylogenetic_statistics(x, profils = None, path_tree = None, path = None, dl_tree = False):
    post.phylogenetic_statistics(x, profils, path_tree, path, dl_tree)
