import pandas as pd
import numpy as np
from scipy.spatial import distance
from scipy import stats
from sklearn.metrics import mutual_info_score
from scipy.stats import fisher_exact
import warnings
warnings.filterwarnings('ignore')


def jaccard(dfx, dfy = None):
    symetry = False
    if dfy is None:
        symetry = True
        dfy = dfx
    jaccard_distance = np.zeros((len(dfx.index), len(dfy.index)))
    for a,i in enumerate(dfx.index): # Éviter les doublons si symétrique
        query = dfx.loc[i]
        for b,j in enumerate(dfy.index):
            if symetry and b<a:
                continue
            score_temp = distance.jaccard(query, dfy.loc[j])
            jaccard_distance[a, b] = score_temp
            if symetry:
                jaccard_distance[b, a] = score_temp
    index=dfx.index 
    columns=dfy.index
    jaccard_distance = pd.DataFrame(jaccard_distance, index=index, columns=columns)
    jaccard_distance = jaccard_distance.fillna(1)
    return jaccard_distance


def hamming(dfx, dfy = None):
    symetry = False
    if dfy is None:
        symetry = True
        dfy = dfx
    hamming_distance = np.zeros((len(dfx.index), len(dfy.index)))
    for a,i in enumerate(dfx.index): # Éviter les doublons si symétrique
        query = dfx.loc[i]
        for b,j in enumerate(dfy.index):
            if symetry and b<a:
                continue
            score_temp = distance.hamming(query, dfy.loc[j])
            hamming_distance[a,b] = score_temp
            if symetry: 
                hamming_distance[b,a] = score_temp
    index=dfx.index 
    columns=dfy.index
    hamming_distance = pd.DataFrame(hamming_distance, index=index, columns=columns)
    hamming_distance = hamming_distance.fillna(1)
    return hamming_distance


def pearson(dfx, dfy = None):
    symetry = False
    if dfy is None:
        dfy = dfx
        symetry=True
    pearson_results = np.zeros((len(dfx.index), len(dfy.index)))
    for a,i in enumerate(dfx.index): # Éviter les doublons si symétrique
        query = dfx.loc[i]
        for b,j in enumerate(dfy.index):
            if symetry and b<a:
                continue
            pearson_correlation = stats.pearsonr(query, dfy.loc[j])
            if not np.isnan(pearson_correlation[0]):
                score_temp = pearson_correlation[0]
                pearson_results[a, b] = score_temp
                if symetry:
                    pearson_results[b, a] = score_temp
            else:
                pearson_results[a, b] = np.nan
                if symetry:
                    pearson_results[b, a] = np.nan
    index=dfx.index 
    columns=dfy.index
    pearson_results = pd.DataFrame(pearson_results, index=index, columns=columns)
    pearson_results = pearson_results.fillna(0)
    return pearson_results


def mi(dfx, dfy = None):
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


def cotransition(tvx, tvy = None, consecutive = True):
    symetry = False
    if tvy is None:
        symetry = True
        tvy = tvx
    cotransition_scores = np.zeros((len(tvx.index), len(tvy.index)))
    p_values = np.zeros((len(tvx.index), len(tvy.index)))
    for a,i in enumerate(tvx.index):
        for b,j in enumerate(tvy.index):
            if symetry and b<a:
                continue
            t1 = 0
            t2 = 0
            c = 0
            d = 0
            k = 0
            v1 = np.array(tvx.loc[i])
            v2 = np.array(tvy.loc[j])
            if consecutive == True :
                t1 = np.count_nonzero(v1)
                t2 = np.count_nonzero(v2)
                nonz = (v1 != 0) & (v2 != 0)
                c = np.sum(v1[nonz]==v2[nonz])
                d = np.count_nonzero(v1[nonz]-v2[nonz])
                k = c - d
            elif consecutive ==  False:
                v1d = np.insert(v1, 0, [0])
                v1dd = np.insert(v1d, 0, [0])
                v2d = np.insert(v2, 0, [0])
                v2dd = np.insert(v2d, 0, [0])
                v1 = np.insert(v1, len(v1), [0, 0])
                v1d = np.insert(v1d, len(v1d), [0])
                v2 = np.insert(v2, len(v2), [0, 0])
                v2d = np.insert(v2d, len(v2d), [0])
                t1 = np.count_nonzero(v1) - np.sum((v1 != 0) & (v1d != 0)) + np.sum((v1 !=0) & (v1d != 0) & (v1dd != 0))
                t2 = np.count_nonzero(v2) - np.sum((v2 != 0) & (v2d != 0)) + np.sum((v2 !=0) & (v2d != 0) & (v2dd != 0))
                nonza = (v1 != 0) & (v2 != 0) & (v1d  == 0) & (v2d == 0)
                nonzb1 = (v1 != 0) & (v2 != 0) & (v1d  == 0)
                nonzb2 = (v1 != 0) & (v2 != 0) & (v2d  == 0)
                nonzc1 = ((v1 != 0) & (v2 != 0)) & ((v1d  != 0) & (v1dd  != 0))
                nonzc2 = ((v1 != 0) & (v2 != 0)) & ((v2d  != 0) & (v2dd  != 0))
                nonz = nonza | (nonzc1 & nonzc2) | (nonzb1 & nonzb2) | (nonzb1 & nonzc2) | (nonzb2 & nonzc1)
                c = np.sum(v1[nonz]==v2[nonz])
                d = np.count_nonzero(v1[nonz]-v2[nonz])
                k = c - d
            if t1 == 0 and t2 == 0:
                cotransition_scores[a, b] = None
                if symetry:
                    cotransition_scores[b, a] = None
            else:
                score_temp = k / (t1 + t2 - abs(k))
                cotransition_scores[a, b] = score_temp
                if symetry:
                    cotransition_scores[b, a] = score_temp

                #tableau de contingence:
            contingency_table = [[abs(k),t1-abs(k)], [t2-abs(k),(len(tvx.columns))-t1-t2+abs(k)]]
            score = fisher_exact(contingency_table, alternative="greater")
            p_values[a, b] = score.pvalue
            if symetry:
                p_values[b, a] = score.pvalue
    cotransition_scores = pd.DataFrame(cotransition_scores, index=tvx.index, columns=tvy.index)
    cotransition_scores = cotransition_scores.fillna(0)
    p_values = pd.DataFrame(p_values, index=tvx.index, columns=tvy.index)
    return cotransition_scores, p_values


def pcs(tvx, tvy, confidence=1.5, penalty=0.6):
    symetry = False
    if tvy is None:
        symetry = True
        tvy = tvx
    pcs_scores = np.zeros((len(tvx.index), len(tvy.index)))
    for a,i in enumerate(tvx.index):
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
            pcs_scores[a, b] = score_temp
            if symetry:
                pcs_scores[b, a] = score_temp
    pcs_scores = pd.DataFrame(pcs_scores, index=tvx.index, columns=tvy.index)
    pcs_scores = pcs_scores.fillna(0)
    return pcs_scores


def SVD_phy(dfx, truncation = 0.5): 
    U, S, V = np.linalg.svd(dfx,False) #SVD de la matrice dfx
    k = int(truncation * len(U)) #nb de colonnes à garder dans U
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
    SVDphy_distance = SVDphy_distance.fillna(1)
    return SVDphy_distance