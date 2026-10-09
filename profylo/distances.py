import pandas as pd
import numpy as np
from scipy.spatial import distance
from scipy import stats
from sklearn.metrics import mutual_info_score
from scipy.stats import fisher_exact
import warnings
from joblib import delayed
from joblib import Parallel
warnings.filterwarnings('ignore')


def jaccard_loop(dfx, dfy, a, i, symetry):
    query = dfx.loc[i]
    row = np.zeros(len(dfy.index))
    for b,j in enumerate(dfy.index):
        if symetry and b<a:
            continue
        row[b] = distance.jaccard(query, dfy.loc[j])
    return a, row

def jaccard(n_job, dfx, dfy = None):
    symetry = False
    if dfy is None:
        symetry = True
        dfy = dfx
    task = [delayed(jaccard_loop)(dfx, dfy, a, i, symetry) for a, i in enumerate(dfx.index)]
    out = Parallel(n_job)(task)
    jaccard_distance = np.zeros((len(dfx.index), len(dfy.index)))
    for a, row in out:
        for b in range(a, len(dfx.index)):
            jaccard_distance[a, b] = row[b]
            if symetry:
                jaccard_distance[b, a] = row[b]
    df_jaccard = pd.DataFrame(jaccard_distance, index=dfx.index , columns=dfy.index).fillna(1)
    return df_jaccard


def hamming_loop(dfx, dfy, a, i, symetry):
    query = dfx.loc[i]
    row = np.zeros(len(dfy.index))
    for b,j in enumerate(dfy.index):
        if symetry and b<a:
            continue
        row[b] = distance.hamming(query, dfy.loc[j])
    return a, row

def hamming(n_job, dfx, dfy = None):
    symetry = False
    if dfy is None:
        symetry = True
        dfy = dfx
    task = [delayed(hamming_loop)(dfx, dfy, a, i, symetry) for a, i in enumerate(dfx.index)]
    out = Parallel(n_job)(task)
    hamming_distance = np.zeros((len(dfx.index), len(dfy.index)))
    for a, row in out:
        for b in range(a, len(dfx.index)):
            hamming_distance[a, b] = row[b]
            if symetry:
                hamming_distance[b, a] = row[b]
    df_hamming = pd.DataFrame(hamming_distance, index=dfx.index , columns=dfy.index).fillna(1)
    return df_hamming


def pearson_loop(dfx, dfy, a, i, symetry):
    query = dfx.loc[i]
    row = np.zeros(len(dfy.index))
    for b,j in enumerate(dfy.index):
        if symetry and b<a:
            continue
        pearson_cor = stats.pearsonr(query, dfy.loc[j])
        if not np.isnan(pearson_cor[0]):
            row[b] = pearson_cor[0]
        else:
            row[b] = np.nan
    return a, row

def pearson(n_job, dfx, dfy = None):
    symetry = False
    if dfy is None:
        symetry = True
        dfy = dfx
    task = [delayed(pearson_loop)(dfx, dfy, a, i, symetry) for a, i in enumerate(dfx.index)]
    out = Parallel(n_job)(task)
    pearson_distance = np.zeros((len(dfx.index), len(dfy.index)))
    for a, row in out:
        for b in range(a, len(dfx.index)):
            pearson_distance[a, b] = row[b]
            if symetry:
                pearson_distance[b, a] = row[b]
    df_jaccard = pd.DataFrame(pearson_distance, index=dfx.index , columns=dfy.index).fillna(0)
    return df_jaccard


def mi_ori(dfx, dfy = None):
    symetry = False
    if dfy is None:
        dfy = dfx
        symetry=True
    mi_distance = np.zeros((len(dfx.index), len(dfy.index)))
    for a,i in enumerate(dfx.index):
        query = dfx.loc[i]
        for b,j in enumerate(dfy.index):
            if symetry and b<a:
                continue
            score_temp = mutual_info_score(query, dfy.loc[j])
            mi_distance[a, b] = score_temp
            if symetry:
                mi_distance[b, a] = score_temp
    index=dfx.index 
    columns=dfy.index
    mi_distance = pd.DataFrame(mi_distance, index=index, columns=columns)
    mi_distance = mi_distance.fillna(0)
    return mi_distance

def mi(dfx, dfy = None):
    npx = dfx.to_numpy()
    symetry = False
    if dfy is None:
        dfy = dfx
        npy = npx
        symetry=True
    mi_distance = np.zeros((len(dfx.index), len(dfy.index)))
    for i in range(len(dfx)):
        query = npx[i]

        p1 = np.sum(query)/len(query)
        if p1 == 1 or p1 ==0:
            entropy_i=0
        else:
            entropy_i = -(p1*np.log(p1) + (1-p1)*np.log(1-p1))

        if symetry:
            loop_range = range(i,len(dfy))
        else:
            loop_range = range(len(dfy))

        for j in loop_range:
            p1 = np.sum(npy[j])/len(npy[j])
            if p1 == 1 or p1 ==0:
                entropy_j = 0
            else:
                entropy_j = -(p1*np.log(p1) + (1-p1)*np.log(1-p1) )
            vsize = len(query)
            p11 = np.sum(np.logical_and(query,npy[j] ))/vsize
            diff = np.logical_xor(query,npy[j])
            p01 = np.sum(np.logical_and(query, diff))/vsize
            p10=  np.sum(np.logical_and(npy[j], diff))/vsize

            p00 = np.sum(np.logical_and(np.logical_not(query),np.logical_not(npy[j] )))/vsize
            cross_ent = -np.sum([x*np.log(x) for x in [p11,p01,p10,p00] if x!=0])
            mi = entropy_i+entropy_j - cross_ent 
            mi = np.round(mi,15)

            mi_distance[i, j] = mi
            if symetry:
                mi_distance[j, i] = mi
    index=dfx.index 
    columns=dfy.index
    mi_distance = pd.DataFrame(mi_distance, index=index, columns=columns)
    mi_distance = mi_distance.fillna(0)
    return mi_distance

def preprocess_transitions(X, consecutive=True):
    if consecutive:
        return X

    nonzero = X != 0

    previous = np.zeros_like(nonzero)
    previous[:, 1:] = nonzero[:, :-1]

    previous2 = np.zeros_like(nonzero)
    previous2[:, 2:] = nonzero[:, :-2]

    second_in_run = nonzero & previous & ~previous2

    return np.where(second_in_run, 0, X)

def concordance(X, Y=None, block_size=1024):
    if Y is None:
        Y = X

    X = np.asarray(X, dtype=np.float32)
    Y = np.asarray(Y, dtype=np.float32)

    n, m = X.shape[0], Y.shape[0]

    K = np.empty((n, m), dtype=np.float32)

    for i in range(0, n, block_size):
        end = min(i + block_size, n)

        K[i:end] = X[i:end] @ Y.T

    return K

def cotransition(tvx, tvy=None, consecutive=True):
    if tvy is None:
        tvy = tvx

    X = tvx.to_numpy(dtype=np.float64)
    Y = tvy.to_numpy(dtype=np.float64)

    X = preprocess_transitions(X, consecutive)
    Y = preprocess_transitions(Y, consecutive)

    print("Pre-processing done...")

    tx = np.count_nonzero(X, axis=1)
    ty = np.count_nonzero(Y, axis=1)

    K = concordance(X, Y)

    denominator = tx[:, None] + ty[None, :] - np.abs(K)

    scores = np.divide(
        K,
        denominator,
        out=np.zeros_like(K),
        where=denominator != 0
    )

    return pd.DataFrame(
        scores,
        index=tvx.index,
        columns=tvy.index
    )

def pcs_loop(tvx, tvy, a, i, symetry, confidence, penalty):
    query = tvx.loc[i]
    row = np.zeros(len(tvy.index))
    for b,j in enumerate(tvy.index):
        if symetry and b<a:
            continue
        match_1 = 0
        mismatch_1 = 0
        match_2 = 0
        mismatch_2 = 0
        tv1 = np.array(tvx.loc[i])
        tv2 = np.array(tvy.loc[j])
        tv1a = np.insert(tv1[:-1], 0, 2) #acceder à la valeur précédente
        tv1b = np.insert(tv1[1:], len(tv1) -1 , 2)  #acceder à la valeur suivante
        tv2a = np.insert(tv2[:-1], 0, 2) 
        tv2b = np.insert(tv2[1:], len(tv2) -1 , 2)
        nonz1 = (tv1 != tv2) & (tv1 != 0)  #True si transition(s) non partagée(s) chez 1
        nonz2 = (tv1 != tv2) & (tv2 != 0)  #True si transition(s) non partagée(s) chez 2
        nonz12= nonz1 & nonz2              #True si transitions différentes chez les deux
        nonzp = (tv1 == tv2) & (tv1 != 0)  #True si transitions partagées
        double_match = (tv1a[nonzp] == 0) & (tv1b[nonzp] == 0) & (tv2a[nonzp] == 0) & (tv2b[nonzp] == 0)
        match_2 = sum(double_match)
        match_1 = sum(nonzp) - match_2
        double_mismatch = sum((tv1a[nonz1] == 0) & (tv1b[nonz1] == 0)) + sum((tv2a[nonz2] == 0) & (tv2b[nonz2] == 0)) - sum((tv1a[nonz12] == 0) & (tv1b[nonz12] == 0) & (tv2a[nonz12] == 0) & (tv2b[nonz12] == 0))
        mismatch_2 = double_mismatch
        mismatch_1 = sum(nonz1) + sum(nonz2) - sum(nonz12) - mismatch_2
        score_temp = (
            (match_1) + (match_2 * confidence) - penalty * (mismatch_1 + mismatch_2 * confidence)
        )
        row[b] = score_temp
    return a, row

def pcs(n_job, tvx, tvy, confidence=1.5, penalty=0.6):
    symetry = False
    if tvy is None:
        symetry = True
        tvy = tvx
    task = [delayed(pcs_loop)(tvx, tvy, a, i, symetry, confidence, penalty) for a, i in enumerate(tvx.index)]
    out = Parallel(n_job)(task)
    pcs_score = np.zeros((len(tvx.index), len(tvy.index)))
    for a, row in out:
        for b in range(a, len(tvx.index)):
            pcs_score[a, b] = row[b]
            if symetry:
                pcs_score[b, a] = row[b]
    df_pcs = pd.DataFrame(pcs_score, index=tvx.index , columns=tvy.index).fillna(0)
    return df_pcs


def SVD_phy(dfx, truncation = 0.5): 
    U, S, V = np.linalg.svd(dfx,False) #SVD de la matrice dfx
    k = int(truncation * np.shape(U)[1]) #nb de colonnes à garder dans U
    U_truncated = U[:, :k] #ajout des colonnes U
    norms = np.linalg.norm(U_truncated, axis=1, keepdims=True)
    U_truncated = U_truncated / norms
    subset_U = U_truncated
    row_labels = dfx.index
    col_labels = dfx.index
    SVDphy_distance = np.zeros((len(subset_U), len(U_truncated)))
    for i in range(0, len(subset_U)):
        for j in range(0, len(U_truncated)):
            SVDphy_distance[i][j] = np.sqrt(np.sum((subset_U[i] - U_truncated[j])**2))
    SVDphy_distance = pd.DataFrame(SVDphy_distance, index=row_labels, columns=col_labels)
    SVDphy_distance = SVDphy_distance.fillna(np.nanmax(SVDphy_distance.to_numpy()))
    return SVDphy_distance
