#!/usr/bin/env python
# coding: utf-8

# In[2]:


import tarfile
import urllib

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, MinMaxScaler, MaxAbsScaler, StandardScaler

from mpl_toolkits.axes_grid1 import make_axes_locatable
# from nuc_tool import CalcNucPositions
from pathlib import Path
from scipy.ndimage import median_filter
import h5py, re
from scipy import stats
import re
import tqdm


# In[3]:


get_ipython().system('pip install openpyxl')


# In[4]:


## standard functions

## function to read fasta file
def read_fasta(filename, name=None):
    with open(filename, 'r') as f:
        content = f.readlines()
        for line in range(len(content)):
            if content[line][0] == '>':
                seq = content[line + 1]
    seq = ''.join([s for s in seq.upper() if s in ['A', 'C', 'G', 'T']])
    return seq

## function to sort LLR traces on summed LLR
def sort_traces(traces, strands, method='sum'):
    from sklearn import decomposition
    new_traces = np.empty_like(traces)
    for strand in np.unique(strands):
        selection = traces[strands == strand]
        if method == 'pca':
            pca = decomposition.PCA(n_components=1)
            pca.fit(traces[strands == strand])
            data = pca.transform(selection)
            i = np.argsort(data[:, 0])
        else:
            i = np.argsort(np.sum(selection, axis=1))
        new_traces[strands == strand] = selection[i]
    return new_traces

## function to find position of motif; a region of interest defined as for example 'acgtcaa'
def find_motif(seq, motif, context=None, format='string'):
    seq = seq.lower()
    motif = motif.lower()
    if context is not None:
        context = context.lower()
    found_bps = seq.replace(motif, motif.upper())
    if context is not None:
        context_bps = seq.replace(context, context.upper())
        found_bps = ''.join([a if np.char.isupper(b) else b for a, b in zip(context_bps, found_bps)])
    if format == 'string':
        return found_bps
    if format == 'numeric':
        return np.asarray([1.0 if bp == bp.upper() else 0 for bp in found_bps])
    if format == 'index':
        return np.asarray([m.start() for m in re.finditer(motif, seq)])

## function to find position of bases that could be methylated
def find_index_methylation(filename, fasta):
    seq = read_fasta(fasta).upper()

    if '5mC' in filename:
        motifsplus = ['Cg', 'gC']
        motifsplus2 = ['gatCg', 'Cgatc', 'gCcagg', 'ccaggC', 'cctggC', 'gCctgg']
        motifsmin = ['Gc', 'cG']
        motifsmin2 = ['ctaGc', 'Gctag', 'cGgtcc', 'ggtccG', 'ggaccG', 'cGgacc']
    elif '6mA' in filename:
        motifsplus = ['A']
        motifsmin = ['T']
    else:
        motifs = None

    new_indexplus = []
    new_indexmin = []

    for motif in motifsplus:
        index = find_motif(seq, motif, format = 'index')
        offset = [i for i, c in enumerate(motif) if c.isupper()]
        new_indexplus.extend(index+offset[0])
    
    for motif in motifsplus2:
        index2 = find_motif(seq, motif, format = 'index')
        offset2 = [i for i, c in enumerate(motif) if c.isupper()]
        index2 = index2 + offset2
        for i in index2:
            new_indexplus.remove(i)
    new_indexplus.sort()

    for motif in motifsmin:
        index = find_motif(seq, motif, format = 'index')
        offset = [i for i, c in enumerate(motif) if c.isupper()]
        new_indexmin.extend(index+offset[0])
        
    for motif in motifsmin2:
        index2 = find_motif(seq, motif, format = 'index')
        offset2 = [i for i, c in enumerate(motif) if c.isupper()]
        index2 = index2 + offset2
        for i in index2:
            new_indexmin.remove(i)
    new_indexmin.sort()

    return np.asarray(new_indexplus), np.asarray(new_indexmin)

## function to read files containing LLR, now only read + or - strand, not both at the same time
def read_files(filename, fasta, strand=''):
    seq = read_fasta(fasta).upper()

    askstrand = strand
    new_indexplus, new_indexmin = find_index_methylation(filename, fasta)

    if '5mC' in filename:
        if askstrand == '+':
            motifs = ['Cg', 'gC']
            motifs2 = ['gatCg', 'Cgatc', 'gCcagg', 'ccaggC', 'cctggC', 'gCctgg']
        elif askstrand == '-':
            motifs = ['Gc', 'cG']
            motifs2 = ['ctaGc', 'Gctag', 'cGgtcc', 'ggtccG', 'ggaccG', 'cGgacc']
        else:
            motifs = ['Cg', 'Gc', 'cG', 'gC']
    elif '6mA' in filename:
        if askstrand == '+':
            motifs = ['A']
        elif askstrand == '-':
            motifs = ['T']
        else:
            motifs = ['A', 'T']
    else:
        motifs = None

    llr, strands = read_stats(filename, seq, length = 4000, motifs = motifs, motifs2 = motifs2, askstrand = askstrand)

    if askstrand == '+':
        new_llr = llr[:, new_indexplus]
    elif askstrand == '-':
        new_llr = llr[:, new_indexmin]

    return llr, new_llr, strands

## function to find llr of bases that could be methylated, index = find_index_methylation (+ or - strand)
def llr_in_range (new_index, fasta, motif, range = None):
    seq = read_fasta(fasta).upper()

    index = find_motif(seq, motif, format = 'index')[0]
    seq = find_motif(seq, motif, format = 'string')
    xlabels = [s + f'\n:\n{i}' if i % 10 == 0 else s for i, s in enumerate(seq)]

    plotrange = [index, index+len(motif)]
    
    if range == None:
        withinrange = new_index[(new_index >= plotrange[0])& (new_index <= plotrange[1])]
    else:
        withinrange = new_index[(new_index >= range[0]) & (new_index <= range[1])]
    
    return plotrange, withinrange


# In[5]:


## To plot normalized sum of LLR (comparing GAL, RAF and no WCE)
# LLRfilename1 = file with LLR of barcode01 (+ or - strand)
# LLRfilename2 = file with LLR of barcode02 (+ or - strand)
# LLRfilename3 = file with LLR of barcode03 (+ or - strand)
# fasta = fasta file of GAL locus
# motif = region of GAL locus you want to look at, defined as for example seqProm
# new_index = array with positions of bases that could have been methylated, defined differently for + or - strand, use find_methylation function to find these positions

def SumPerSiteNewVersion(LLRfilename):
    Sumllr = np.sum(-LLRfilename, axis = 0)
    scaler = MaxAbsScaler()
    Sumllr = scaler.fit_transform([[x] for x in Sumllr])

    return Sumllr

def plotSumPerSite(LLRfilename, fasta, new_index, motif):
    plotrange, withinrange = llr_in_range(new_index, fasta, motif)
    
    Sumllr = SumPerSiteNewVersion(LLRfilename[:,withinrange])
    
    plt.scatter(withinrange, Sumllr)
    
def plotSumPerSiteComparingConditions(LLRfilename1, LLRfilename2, LLRfilename3, fasta, new_index, motif, outputdir, namefig):
    plotrange, withinrange = llr_in_range(new_index, fasta, motif)
    
    Sumllr1 = SumPerSiteNewVersion(LLRfilename1[:,withinrange])
    Sumllr2 = SumPerSiteNewVersion(LLRfilename2[:,withinrange])
    Sumllr3 = SumPerSiteNewVersion(LLRfilename3[:,withinrange])

    Sum1min2 = Sumllr1-Sumllr2
    Sum1min3 = Sumllr1-Sumllr3

    fig1 = plt.figure(figsize = (14, 10), dpi = 200)

    plt.step(withinrange, Sumllr1, label = 'GAL', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
    plt.step(withinrange, Sumllr2, label = 'RAF', linewidth = 1, linestyle = '--') #, s = 50, marker = 's')
    plt.step(withinrange, Sumllr3, label = 'No WCE', linewidth = 1, linestyle = '-.') #, s = 50, marker = 'P')
    
    plt.ylabel('average LLR')
    plt.ylim((-1,1))
    plt.xlabel('i')
    plt.legend()
    plt.grid(True)
    
    plt.savefig(outputdir + namefig, dpi=1200, transparent=True)
    plt.show()


# In[6]:


## K-means clustering of LLR of sequenced DNA
# LLRfilename = file with LLR of any barcode (+ or - strand)
# fasta = fasta file of GAL locus
# motif = region of GAL locus you want to look at, defined as for example seqProm
# new_index = array with positions of bases that could have been methylated, defined differently for + or - strand, use find_methylation function to find these positions

def clustering_kmeans_newversion (subllr, n_clusters):   # k-means  
    preprocessor = Pipeline (
        [
            ('scaler', MaxAbsScaler()),
            ('pca', PCA(n_components = 2, random_state = 42)),
        ]
    )

    clusterer = Pipeline(
        [
            (
                "kmeans",
                KMeans(
                    n_clusters = n_clusters,
                    init ='k-means++',
                    n_init = 50,
                    max_iter = 500,
                    random_state = 42,
                    ),
                ),
        ]
    )

    pipe = Pipeline(
        [
            ("preprocessor", preprocessor),
            ("clusterer", clusterer)
        ]
    )

    pipe.fit(subllr)
    preprocessed_data = pipe["preprocessor"].transform(subllr)
    predicted_labels = pipe["clusterer"]["kmeans"].labels_

    pcadf = pd.DataFrame(
        pipe["preprocessor"].transform(subllr),
        columns = ["component_1", "component_2"],
    )
    pcadf["predicted_cluster"] = pipe["clusterer"]["kmeans"].labels_

    return preprocessed_data, predicted_labels

def clustering_kmeans_determinantNewVersion(LLRfilename, fasta, motif, new_index, outputdir, namefig):  
    plotrange, withinrange = llr_in_range (new_index, fasta, motif, range = None)
    
    subllr = LLRfilename[:, withinrange]

    silhouette_scores = -100
    n_clusterfinal = 0
    
    for i in range(2,11):
        n_clusters=int(i)

        preprocessed_data, predicted_labels = clustering_kmeans_newversion (subllr, n_clusters)
        silhouette_coef = silhouette_score(preprocessed_data, predicted_labels)

        if silhouette_coef > silhouette_scores and silhouette_coef <= 1:
            silhouette_scores = silhouette_coef
            n_clusterfinal = i

    preprocessor = Pipeline (
        [
            ('scaler', MaxAbsScaler()),
            ('pca', PCA(n_components = 2, random_state = 42)),
        ]
    )

    clusterer = Pipeline(
        [
            (
                "kmeans",
                KMeans(
                    n_clusters=n_clusterfinal,
                    init ='k-means++',
                    n_init = 50,
                    max_iter = 500,
                    random_state = 42,
                    ),
                ),
        ]
    )

    pipe = Pipeline(
        [
            ("preprocessor", preprocessor),
            ("clusterer", clusterer)
        ]
    )

    pipe.fit(subllr)
    preprocessed_data = pipe["preprocessor"].transform(subllr)
    predicted_labels = pipe["clusterer"]["kmeans"].labels_
    

    pcadf = pd.DataFrame(
        pipe["preprocessor"].transform(subllr),
        columns = ["component_1", "component_2"],
    )
    #pcadf["predicted_cluster"] = pipe["clusterer"]["kmeans"].labels_ 
    
    plt.figure(figsize = (8,8)) 
    #plt.style.use("fivethirtyeight")
    
    scat = sns.scatterplot(
        x = "component_1",
        y = "component_2",
        s = 50,
        data = pcadf,
        hue = predicted_labels,
        palette = "Set2",
    )

    plt.legend(bbox_to_anchor = (1.05,1), loc = 2, borderaxespad = 0.0)
    plt.legend()
    
    plt.savefig(outputdir + namefig, dpi=1200, transparent=True)
    plt.show()
    
    return n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange


# In[7]:


## Hierarchical clustering: AgglomerativeClustering
# LLRfilename = file with LLR of any barcode (+ or - strand)
# fasta = fasta file of GAL locus
# motif = region of GAL locus you want to look at, defined as for example seqProm
# new_index = array with positions of bases that could have been methylated, defined differently for + or - strand, use find_methylation function to find these positions

from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as shc


def Agglomerative_clustering(pcadf, n_cluster): 
    hierarchical_cluster = AgglomerativeClustering(n_clusters = n_cluster, affinity = 'euclidean', linkage = 'ward')
    predicted_labels = hierarchical_cluster.fit_predict(pcadf)
    
    return predicted_labels
    
def clustering_agglomerative_determinant(LLRfilename, fasta, motif, new_index, outputdir, namefig):
    plotrange, withinrange = llr_in_range(new_index, fasta, motif)
    
    subllr = LLRfilename[:, withinrange]

    preprocessor = Pipeline (
        [
            ('scaler', MaxAbsScaler()),
            ('pca', PCA(n_components = 2, random_state = 42)),
        ]
    )

    pipe = Pipeline(
        [
            ("preprocessor", preprocessor),
        ]
    )

    pipe.fit(subllr)
    preprocessed_data = pipe["preprocessor"].transform(subllr)

    pcadf = pd.DataFrame(
        pipe["preprocessor"].transform(subllr),
        columns = ["component_1", "component_2"],
    )

    silhouette_scores = -100
    n_clusterfinal = 0
    
    for i in range(2,11):
        n_cluster = int(i)

        predicted_labels = Agglomerative_clustering(pcadf, n_cluster)
        silhouette_coef = silhouette_score(pcadf, predicted_labels)

        if silhouette_coef > silhouette_scores and silhouette_coef <= 1:
            silhouette_scores = silhouette_coef
            n_clusterfinal = i

    predicted_labels = Agglomerative_clustering(pcadf, n_clusterfinal)
    
    #plt.style.use("fivethirtyeight")
    plt.figure(figsize = (8,8))

    scat = sns.scatterplot(
        x = "component_1",
        y = "component_2",
        s = 50,
        data = pcadf,
        hue = predicted_labels,
        palette = "Set2",
    )

    plt.legend(bbox_to_anchor = (1.05,1), loc = 2, borderaxespad = 0.0)
    plt.legend()
    
    plt.savefig(outputdir + namefig, dpi=1200, transparent=True)
    plt.show()
    
    return n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange

def plot_dendrogram(LLRfilename, pcadf, outputdir, namefig):
    plt.style.use('default')
    
    Z = shc.linkage(pcadf, method ='ward', metric='euclidean')
    
    plt.figure(figsize=(14,10), dpi = 200)
    info = shc.dendrogram(Z, leaf_font_size=8, orientation='left')
    #plt.xlim((150,0))
    plt.savefig(outputdir + namefig, dpi=1200, transparent=True)
    plt.show()
    
    new_direction = []
    for i in info['leaves']:
        new_direction.append(LLRfilename[i])
    new_direction = np.asarray(new_direction)
    
    return new_direction
    


# In[8]:


## Plotting clusters and LLR sorted on cluster + function to split clusters into a dictionary
# LLRfilename = file with LLR of any barcode (+ or - strand)
# fasta = fasta file of GAL locus
# motif = region of GAL locus you want to look at, defined as for example seqProm
# index = array with positions of bases that could have been methylated, defined differently for + or - strand, use find_methylation function to find these positions
# clusteringmethod = 'kmeans' or 'agglomerative'

def sort_clustersNewVersion(ncluster, predicted_labels, LLRfilename):
    new_clusters = []
    
    for cluster in range(0,ncluster):
        selection = np.where(predicted_labels == cluster)
        new_clusters.extend(LLRfilename[selection])
    new_clusters = np.asarray(new_clusters)
    return new_clusters

def split_clusters(LLRfilename, ncluster, predicted_labels):
    dictCluster = {}

    for cluster in range(0,ncluster):
        selection = np.where(predicted_labels == cluster)
        dictCluster["Cluster{0}".format(cluster)] = LLRfilename[selection]

    return dictCluster

def plot_clustered_llrNewVersion(LLRfilename, fasta, motif, index, outputdir, namefig, clusteringmethod = ''):
    if clusteringmethod == 'kmeans':
        ncluster,score, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(LLRfilename, fasta, motif, index)
    elif clusteringmethod == 'agglomerative':
        ncluster,score, pcadf, predicted_labels, plotrange = clustering_agglomerative_determinant(LLRfilename, fasta, motif, index)
    
    sorted_clusters = sort_clustersNewVersion(ncluster, predicted_labels, LLRfilename)

    colorrange = 1.0
    title = 'Test' # still have to change
    plt.figure(figsize=(8,8),dpi=200)
    plt.text(0, 1.0, title, horizontalalignment = 'left', verticalalignment = 'bottom', transform = plt.gca().transAxes)

    ax = plt.gca()

    im = ax.imshow(-sorted_clusters, cmap='bwr', vmin = -colorrange, vmax = colorrange, origin = 'lower')
    plt.ylabel(f'molecule#')
    plt.xlim((plotrange[0], plotrange[1]))
    plt.xlabel('i')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size = "1%", pad=0.1)
    cbar = plt.colorbar(im, cax = cax)
    cbar.set_label('LLR', rotation = 270)
    
    plt.savefig(outputdir + namefig, dpi=1200, transparent=True)
    plt.show()
    
def plot_llr(LLRfilename, strands, fasta, motif, new_index, outputdir, namefig, method = '', predicted_labels = None, ncluster = None):
    plt.style.use('default')
    
    plotrange, withinrange = llr_in_range (new_index, fasta, motif, range = None)
    
    if method == 'sum':
        sorted_llr = sort_traces(LLRfilename, strands)
    elif method == 'clusters':
        sorted_llr = sort_clustersNewVersion(ncluster, predicted_labels, LLRfilename)
    elif method == 'dendrogram':
        sorted_llr = LLRfilename # input LLRfilename is new_direction of plot_dendrogram (block 6)
    else: 
        sorted_llr = LLRfilename

    colorrange = 1.0
    title = '' # still have to change
    plt.figure(figsize = (14,10) ,dpi=200) # for whole GAL locus 14,10
    plt.text(0, 1.0, title, horizontalalignment = 'left', verticalalignment = 'bottom', transform = plt.gca().transAxes)

    ax = plt.gca()
  
    im = ax.imshow(-sorted_llr, cmap='bwr', vmin = -colorrange, vmax = colorrange, origin = 'lower')
    plt.ylabel(f'molecule#')
    plt.xlim((plotrange[0], plotrange[1]))
    plt.xlabel('i')
    #plt.grid(False)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size = "1%", pad=0.1) # size = 1% and pad = 0.1
    cbar = plt.colorbar(im, cax = cax)
    cbar.set_label('LLR', rotation = 270)
    
    plt.savefig(outputdir + namefig, dpi=1200, transparent=True) #, edgecolor='black') # dpi = 1200
    plt.show()
    


# In[9]:


## plotting voilinplot, only plot + or - strand
# LLRfilename = file with LLR of any barcode (+ or - strand)
# fasta = fasta file of GAL locus
# motif = region of GAL locus you want to look at, defined as for example seqProm
# askstrand = '+' if you look at + strand, '-' if you look at - strand

def create_voilinplot(LLRfilename, fasta, motif, askstrand, outputdir, namefig):
    seq = read_fasta(fasta).upper()
    
    colorrange = 1.0
    plt.figure(figsize=(40, 2), dpi=100)
    #title = filename.split(r'/')[-2] + f'\n\n{filename.split(".")[-3]} {label}'
    plt.text(0, 1.0, 'Test', horizontalalignment='left', verticalalignment='bottom', transform=plt.gca().transAxes)
    ax = plt.gca()
    
    import matplotlib.patches as mpatches

    index = find_motif(seq, motif, format='index')[0]
    seq = find_motif(seq, motif, format='string')
    xlabels = [s + f'\n:\n{i+1}' if i % 10 == 0 else s for i, s in enumerate(seq)]
    #xlabels = [s + f'\n:\n{i}' if i % 10 == 0 else s for i, s in enumerate(seq)]

    plotrange= [index, index+len(motif)]
    #xlabels = xlabels[plotrange]
    xlabels = xlabels[plotrange[0]:plotrange[1]]

    labels = []

    if askstrand == '+':
        s = 0
        color = 'red'
    elif askstrand == '-':
        s = 1
        color = 'blue'

    subllr = LLRfilename[:, plotrange[0]:plotrange[1]]
    violin_parts = plt.violinplot(-subllr, showextrema=False)
    for pc in violin_parts['bodies']:
        pc.set_facecolor(color)
        pc.set_edgecolor('black')
            
    labels.append((mpatches.Patch(color=color), f'strand {s}'))

    plt.legend(*zip(*labels), loc=2)
    plt.xticks(np.arange(len(xlabels)) + 1, xlabels)

    plt.ylim((-10, 10))
    plt.ylabel(f'-LLR')
    plt.tight_layout()
    
    plt.savefig(outputdir + namefig, dpi=1200)
    plt.show()
    


# In[10]:


## motifs

# seq601 = 'CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGCAAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT'
seqGAL = 'TCAGGAAAATCTGTAGACAATCTTGGACCCGTAAGTTTCACCGTTTTTCAAGGTTACACAATCTTTCCAGTTCTCTTGATTGATAGCATCAATGTATCTACCAGGCTCAATTGCAAAACCTTGTCTTGCTTCGTAACCAGCAGACAAGAAATCACCGGTATAAAATTGATAAGTTGGCTCTGTACTTAAAACTTCTAATGTAATATTGGAATCGGGATGAAAAGCCTTGACAATAAGCGTCAATTCATTGTTTAGAGTATTGATTTGACTTGGCTTAGCATTTTCATCCACCACAAAACAACAATCAAACTGGGGATTTTTGGGGCCTAAGACCGTTGGCTTTGTAGAGTTAAAGGTAGCAATTTCTCTATCGACGATATTACCCGTAGGAATCATGTTTTTGTCGACATCAACAGATTTTTTTGAACGCACCATAATCTCCGTACCCTCAATAGTGTCTCCATATGGCTTGTTCAGATTGAAATAACTATGATTTGTTAAATTTATTGGCGTCGCTTCACCAGCAGTCAATTTACCTTTATATACCATTTCCAAACTTTTTTGGGCAACGTTCACAGTATACTGTATGGTTACCAATAGATCACCTGGAAATTCGGTGTCCTTCTCATTATCTATCAGCATGTACTCGGCGGTAAAAACATCCTTTGAAGGATTTTGAATGATGGGTCCCAAAAATCTTTTTCTGTGGAAAGAACCGATACTACTATGATTCGCATTAACGCCGTTATTAACGGTTAACTGATAGTCTTTGTTGCATAAACTAAACTTACCCTTCGAAATACGATTAGCATACCTGCCGATCGTGGCGCCTATATAAGCACTATCAGGATTCAAATACCCTTCCTCATTTTCATAGCCAAGAACAACTGATTGTCCGTTCACTTTCAGGTCAACAATGCTGGCGCCCAAATTGGCAAACGTGGCTTGAAATCTGGTGCCGGCACCAATAGTCACAAATCTTGCGTCATAACGCATATCTTCAGCGGAAAATCTGGCCTCGACACCCCTTAACTGGTAACCAAAAGGATTCTCAGTAGTCCATTTCCATAAATCCTTGCAGGAGTCTTCAACCTGCAACTCGGTCTGCCATTTCAGTTCGCGTTTGGCCCTATCTGGTTTAGCCGTCAAGTTCAAAACATCACCTGCTCTTCTGCCCGTAACTTTGTATGGAAGATCAATACCAGAAGCTTTGCAGAATGCATGATAAACTTCAAAAACTGTAGAACCTTTACCGGAACCCAAGTTCCACTCACGACACAAACCTTCATTTTCATTGTAGGCCTCTAGGTATTGCAGGGCTGCAATATGACCTTTTGCTAGATCAACTACGTGGATATAATCCCTGATCGGGGTACCATCTCTGGAATCATAATCGTCTCCGAAGATGTAAAGCTTCTCGCGCCTACCAACAGCTACTTGAGCCATATATGGCAACAAATTGTTTGGTATACCTAGCGGATCTTCTCCGATTAATCCAGAGGGATGTGCGCCAATTGGGTTAAAATAACGCAAGATAGCAAACTTCCAACTTTTTTTGTCGCTATTGTAAAGATCATTCAAGATATTCTCAATGGCGTATTTCGTATGACCATACGGATTAGTAGGCCCTAAGGGACATTCTTCTGGGATAGGAATCATATTTGGGAATCTCGTAGCATCACCATAGACAGTAGCAGAAGATGAAAAAACAAATTTGGAAACGTTGTATTGTTGCATTAACTCTAATAAAACGACAGTTCCCAAAATGTTATTGTGATAGTATCTCAGCGGGATTTGTGTAGATTCACCTACAGCCTTTAAACCAGCAAAGTGAATTACCGAATCAATTTTATATTCTTTGAAAACCTTTTCCAGACCTTTTCGGTCACACAAATCAACCTCATAGAAGGGAATGTGATGCTTGGTCAAGACCTCTAACCTGGCTACAGAATCATAAGTTGAATTCGACAGGTTATCAGCAACAACACAGTCATATCCATTCTCAATTAGCTCTACCACAGTGTGTGAACCAATGTATCCAGCACCACCTGTAACCAAAACAATTTTAGAAGTACTTTCACTTTGTAACTGAGCTGTCATTTATATTGAATTTTCAAAAATTCTTACTTTTTTTTTGGATGGACGCAAAGAAGTTTAATAATCATATTACATGGCATTACCACCATATACATATCCATATCTAATCTTACTTATATGTTGTGGAAATGTAAAGAGCCCCATTATCTTAGCCTAAAAAAACCTTCTCTTTGGAACTTTCAGTAATACGCTTAACTGCTCATTGCTATATTGAAGTACGGATTAGAAGCCGCCGAGCGGGCGACAGCCCTCCGACGGAAGACTCTCCTCCGTGCGTCCTCGTCTTCACCGGTCGCGTTCCTGAAACGCAGATGTGCCTCGCGCCGCACTGCTCCGAACAATAAAGATTCTACAATACTAGCTTTTATGGTTATGAAGAGGAAAAATTGGCAGTAACCTGGCCCCACAAACCTTCAAATTAACGAATCAAATTAACAACCATAGGATGATAATGCGATTAGTTTTTTAGCCTTATTTCTGGGGTAATTAATCAGCGAAGCGATGATTTTTGATCTATTAACAGATATATAAATGGAAAAGCTGCATAACCACTTTAACTAATACTTTCAACATTTTCAGTTTGTATTACTTCTTATTCAAATGTCATAAAAGTATCAACAAAAAATTGTTAATATACCTCTATACTTTAACGTCAAGGAGAAAAAACTATAATGACTAAATCTCATTCAGAAGAAGTGATTGTACCTGAGTTCAATTCTAGCGCAAAGGAATTACCAAGACCATTGGCCGAAAAGTGCCCGAGCATAATTAAGAAATTTATAAGCGCTTATGATGCTAAACCGGATTTTGTTGCTAGATCGCCTGGTAGAGTCAATCTAATTGGTGAACATATTGATTATTGTGACTTCTCGGTTTTACCTTTAGCTATTGATTTTGATATGCTTTGCGCCGTCAAAGTTTTGAACGAGAAAAATCCATCCATTACCTTAATAAATGCTGATCCCAAATTTGCTCAAAGGAAGTTCGATTTGCCGTTGGACGGTTCTTATGTCACAATTGATCCTTCTGTGTCGGACTGGTCTAATTACTTTAAATGTGGTCTCCATGTTGCTCACTCTTTTCTAAAGAAACTTGCACCGGAAAGGTTTGCCAGTGCTCCTCTGGCCGGGCTGCAAGTCTTCTGTGAGGGTGATGTACCAACTGGCAGTGGATTGTCTTCTTCGGCCGCATTCATTTGTGCCGTTGCTTTAGCTGTTGTTAAAGCGAATATGGGCCCTGGTTATCATATGTCCAAGCAAAATTTAATGCGTATTACGGTCGTTGCAGAACATTATGTTGGTGTTAACAATGGCGGTATGGATCAGGCTGCCTCTGTTTGCGGTGAGGAAGATCATGCTCTATACGTTGAGTTCAAACCGCAGTTGAAGGCTACTCCGTTTAAATTTCCGCAATTAAAAAACCATGAAATTAGCTTTGTTATTGCGAACACCCTTGTTGTATCTAACAAGTTTGAAACCGCCCCAACCAACTATAATTTAAGAGTGGTAGAAGTCACTACAGCTGCAAATGTTTTAGCTGCCACGTACGGTGTTGTTTTACTTTCTGGAAAAGAAGGATCGAGCACGAATAAAGGTAATCTAAGAGATTTCATGAACGTTTATTATGCCAGATATCACAACATTTCCACACCCTGGAACGGCGATATTGAATCCGGCATCGAACGGTTAACAAAGATGCTAGTACTAGTTGAAGAGTCTCTCGCCAATAAGAAACAGGGCTTTAGTGTTGACGATGTCGCACAATCCTTGAATTGTTCTCGCGAAGAATTCACAAGAGACTACTTAACAACATCTCCAGTGAGATTTCAAGTCTTAAAGCTATATCAGAGGGCTAAGCATGTGTATTCTGAATCTTTAAGAGTCTTGAAGGCTGTGAAATTAATGACTACAGCGAGCTTTACTGCCGACGAAGACTTTTTCAAGCAATTTGGTGCCTTGATGAACGAGTCTCAAGCTTCTTGCGATAAACTTTACGAATGTTCTTGTCCAGAGATTGACAAAATTTGTTCCATTGCTTTGTCAAATGGATCATATGGTTCCCGTTTGACCGGAGCTGGCTGGGGTGGTTGTACTGTTCACTTGGTTCCAGGGGGCCCAAATGGCAACATAGAAAAGGTAAAAGAAGCCCTTGCCAATGAGTTCTACAAGGTCAAGTACCCTAAGATCACTGATGCTGAGCTAGAAAATGCTATCATCGTCTCTAAACCAGCATTGGGCAGCTGTCTATATGAATTATAA'
seqProm = 'TTATATTGAATTTTCAAAAATTCTTACTTTTTTTTTGGATGGACGCAAAGAAGTTTAATAATCATATTACATGGCATTACCACCATATACATATCCATATCTAATCTTACTTATATGTTGTGGAAATGTAAAGAGCCCCATTATCTTAGCCTAAAAAAACCTTCTCTTTGGAACTTTCAGTAATACGCTTAACTGCTCATTGCTATATTGAAGTACGGATTAGAAGCCGCCGAGCGGGCGACAGCCCTCCGACGGAAGACTCTCCTCCGTGCGTCCTCGTCTTCACCGGTCGCGTTCCTGAAACGCAGATGTGCCTCGCGCCGCACTGCTCCGAACAATAAAGATTCTACAATACTAGCTTTTATGGTTATGAAGAGGAAAAATTGGCAGTAACCTGGCCCCACAAACCTTCAAATTAACGAATCAAATTAACAACCATAGGATGATAATGCGATTAGTTTTTTAGCCTTATTTCTGGGGTAATTAATCAGCGAAGCGATGATTTTTGATCTATTAACAGATATATAAATGGAAAAGCTGCATAACCACTTTAACTAATACTTTCAACATTTTCAGTTTGTATTACTTCTTATTCAAATGTCATAAAAGTATCAACAAAAAATTGTTAATATACCTCTATACTTTAACGTCAAGGAGAAAAAACTATA'
seqProm_Gal1 = 'CGGATTAGAAGCCGCCGAGCGGGCGACAGCCCTCCGACGGAAGACTCTCCTCCGTGCGTCCTCGTCTTCACCGGTCGCGTTCCTGAAACGCAGATGTGCCTCGCGCCGCACTGCTCCGAACAATAAAGATTCTACAATACTAGCTTTTATGGTTATGAAGAGGAAAAATTGGCAGTAACCTGGCCCCACAAACCTTCAAATTAACGAATCAAATTAACAACCATAGGATGATAATGCGATTAGTTTTTTAGCCTTATTTCTGGGGTAATTAATCAGCGAAGCGATGATTTTTGATCTATTAACAGATATATAAATGGAAAAGCTGCATAACCACTTTAACTAATACTTTCAACATTTTCAGTTTGTATTACTTCTTATTCAAATGTCATAAAAGTATCAACAAAAAATTGTTAATATACCTCTATACTTTAACGTCAAGGAGAAAAAACTATAATGACTAAATCTCATTCAGAAGAAGTGATTGTACCTGAGTTCAATTCTAGCGCAAAGGAATTACCAAGACCATTGGCCGAAAAGTGCCCGAGCATAATTAAGAAATTTATAAGCGCTTATGATGCTAAACCGGATTTTGTTGCTAGATCGCCTGGTAGAGTCAATCTAATTGGTGAACATATTGATTATTGTGACTTCTCGGTTTTACCTTTAGCTATTGATTTTGATATGCTTTGCGCCGTCAAAGTTTTGAACGAGAAAAATCCATCCATTACCTTAATAAATGCTGATCCCAAATTTGCTCAAAGGAAGTTCGATTTGCCGTTGGACGGTTCTTATGTCACAATTGATCCTTCTGTGTCGGACTGGTCTAATTACTTTAAATGTGGTCTCCATGTTGCTCACTCTTTTCTAAAGAAACTTGCACCGGAAAGGTTTGCCAGTGCTCCTCTGGCCGGGCTGCAAGTCTTCTGTGAGGGTGATGTACCAACTGGCAGTGGATTGTCTTCTTCGGCCGCATTCATTTGTGCCGTTGCTTTAGCTGTTGTTAAAGCGAATATGGGCCCTGGTTATCATATGTCCAAGCAAAATTTAATGCGTATTACGGTCGTTGCAGAACATTATGTTGGTGTTAACAATGGCGGTATGGATCAGGCTGCCTCTGTTTGCGGTGAGGAAGATCATGCTCTATACGTTGAGTTCAAACCGCAGTTGAAGGCTACTCCGTTTAAATTTCCGCAATTAAAAAACCATGAAATTAGCTTTGTTATTGCGAACACCCTTGTTGTATCTAACAAGTTTGAAACCGCCCCAACCAACTATAATTTAAGAGTGGTAGAAGTCACTACAGCTGCAAATGTTTTAGCTGCCACGTACGGTGTTGTTTTACTTTCTGGAAAAGAAGGATCGAGCACGAATAAAGGTAATCTAAGAGATTTCATGAACGTTTATTATGCCAGATATCACAACATTTCCACACCCTGGAACGGCGATATTGAATCCGGCATCGAACGGTTAACAAAGATGCTAGTACTAGTTGAAGAGTCTCTCGCCAATAAGAAACAGGGCTTTAGTGTTGACGATGTCGCACAATCCTTGAATTGTTCTCGCGAAGAATTCACAAGAGACTACTTAACAACATCTCCAGTGAGATTTCAAGTCTTAAAGCTATATCAGAGGGCTAAGCATGTGTATTCTGAATCTTTAAGAGTCTTGAAGGCTGTGAAATTAATGACTACAGCGAGCTTTACTGCCGACGAAGACTTTTTCAAGCAATTTGGTGCCTTGATGAACGAGTCTCAAGCTTCTTGCGATAAACTTTACGAATGTTCTTGTCCAGAGATTGACAAAATTTGTTCCATTGCTTTGTCAAATGGATCATATGGTTCCCGTTTGACCGGAGCTGGCTGGGGTGGTTGTACTGTTCACTTGGTTCCAGGGGGCCCAAATGGCAACATAGAAAAGGTAAAAGAAGCCCTTGCCAATGAGTTCTACAAGGTCAAGTACCCTAAGATCACTGATGCTGAGCTAGAAAATGCTATCATCGTCTCTAAACCAGCATTGGGCAGCTGTCTATATGAATTATAA'
seqProm_Gal10 = 'TCAGGAAAATCTGTAGACAATCTTGGACCCGTAAGTTTCACCGTTTTTCAAGGTTACACAATCTTTCCAGTTCTCTTGATTGATAGCATCAATGTATCTACCAGGCTCAATTGCAAAACCTTGTCTTGCTTCGTAACCAGCAGACAAGAAATCACCGGTATAAAATTGATAAGTTGGCTCTGTACTTAAAACTTCTAATGTAATATTGGAATCGGGATGAAAAGCCTTGACAATAAGCGTCAATTCATTGTTTAGAGTATTGATTTGACTTGGCTTAGCATTTTCATCCACCACAAAACAACAATCAAACTGGGGATTTTTGGGGCCTAAGACCGTTGGCTTTGTAGAGTTAAAGGTAGCAATTTCTCTATCGACGATATTACCCGTAGGAATCATGTTTTTGTCGACATCAACAGATTTTTTTGAACGCACCATAATCTCCGTACCCTCAATAGTGTCTCCATATGGCTTGTTCAGATTGAAATAACTATGATTTGTTAAATTTATTGGCGTCGCTTCACCAGCAGTCAATTTACCTTTATATACCATTTCCAAACTTTTTTGGGCAACGTTCACAGTATACTGTATGGTTACCAATAGATCACCTGGAAATTCGGTGTCCTTCTCATTATCTATCAGCATGTACTCGGCGGTAAAAACATCCTTTGAAGGATTTTGAATGATGGGTCCCAAAAATCTTTTTCTGTGGAAAGAACCGATACTACTATGATTCGCATTAACGCCGTTATTAACGGTTAACTGATAGTCTTTGTTGCATAAACTAAACTTACCCTTCGAAATACGATTAGCATACCTGCCGATCGTGGCGCCTATATAAGCACTATCAGGATTCAAATACCCTTCCTCATTTTCATAGCCAAGAACAACTGATTGTCCGTTCACTTTCAGGTCAACAATGCTGGCGCCCAAATTGGCAAACGTGGCTTGAAATCTGGTGCCGGCACCAATAGTCACAAATCTTGCGTCATAACGCATATCTTCAGCGGAAAATCTGGCCTCGACACCCCTTAACTGGTAACCAAAAGGATTCTCAGTAGTCCATTTCCATAAATCCTTGCAGGAGTCTTCAACCTGCAACTCGGTCTGCCATTTCAGTTCGCGTTTGGCCCTATCTGGTTTAGCCGTCAAGTTCAAAACATCACCTGCTCTTCTGCCCGTAACTTTGTATGGAAGATCAATACCAGAAGCTTTGCAGAATGCATGATAAACTTCAAAAACTGTAGAACCTTTACCGGAACCCAAGTTCCACTCACGACACAAACCTTCATTTTCATTGTAGGCCTCTAGGTATTGCAGGGCTGCAATATGACCTTTTGCTAGATCAACTACGTGGATATAATCCCTGATCGGGGTACCATCTCTGGAATCATAATCGTCTCCGAAGATGTAAAGCTTCTCGCGCCTACCAACAGCTACTTGAGCCATATATGGCAACAAATTGTTTGGTATACCTAGCGGATCTTCTCCGATTAATCCAGAGGGATGTGCGCCAATTGGGTTAAAATAACGCAAGATAGCAAACTTCCAACTTTTTTTGTCGCTATTGTAAAGATCATTCAAGATATTCTCAATGGCGTATTTCGTATGACCATACGGATTAGTAGGCCCTAAGGGACATTCTTCTGGGATAGGAATCATATTTGGGAATCTCGTAGCATCACCATAGACAGTAGCAGAAGATGAAAAAACAAATTTGGAAACGTTGTATTGTTGCATTAACTCTAATAAAACGACAGTTCCCAAAATGTTATTGTGATAGTATCTCAGCGGGATTTGTGTAGATTCACCTACAGCCTTTAAACCAGCAAAGTGAATTACCGAATCAATTTTATATTCTTTGAAAACCTTTTCCAGACCTTTTCGGTCACACAAATCAACCTCATAGAAGGGAATGTGATGCTTGGTCAAGACCTCTAACCTGGCTACAGAATCATAAGTTGAATTCGACAGGTTATCAGCAACAACACAGTCATATCCATTCTCAATTAGCTCTACCACAGTGTGTGAACCAATGTATCCAGCACCACCTGTAACCAAAACAATTTTAGAAGTACTTTCACTTTGTAACTGAGCTGTCATTTATATTGAATTTTCAAAAATTCTTACTTTTTTTTTGGATGGACGCAAAGAAGTTTAATAATCATATTACATGGCATTACCACCATATACATATCCATATCTAATCTTACTTATATGTTGTGGAAATGTAAAGAGCCCCATTATCTTAGCCTAAAAAAACCTTCTCTTTGGAACTTTCAGTAATACGCTTAACTGCTCATTGCTATATTGAAGTACGGATTAGAAGCCGCCGAGCGGGCGACAGCCCTCCGACGGAAGACTCTCCTCCGTGCGTCCTCGTCTTCACCGGTCGCGTTCCTGAAACGCAGATGTGCCTCGCGCCGCACTGCTCCGAACAATAAAGATTCTACAATACTAGCTTTTATGGTTATGAAGAGGAAAAATTGGCAGTAACCTGGCCCCACAAACCTTCAAATTAACGAATCAAATTAACAACCATAGGATGATAATGCGATTAGTTTTTTAGCCTTATTTCTGGGGTAATTAATCAGCGAAGCGATGATTTTTGATCTATTAACAGATATATAAATGGAAAAGCTGCATAACCACTTTAACTAATACTTTCAACATTTTCAGTTTGTATTACTTCTTATTCAAATGTCATAAAAGTATCAACAAAAAATTGTTAATATACCTCTATACTTTAACGTCAAGGAG'
seqGal4_BS1 = 'CGGATTAGAAGCCGCCG'
seqGal4_BS2 = 'CGGGCGACAGCCCTCCG'
seqGal4_BS3 = 'CGGAAGACTCTCCTCCG'
seqGal4_BS4 = 'CTCGCGCCGCACTGCTCCG'
seqGal4_UAS = 'CGGATTAGAAGCCGCCGAGCGGGCGACAGCCCTCCGACGGAAGACTCTCCTCCGTGCGTCCTCGTCTTCACCGGTCGCGTTCCTGAAACGCAGATGTGCCTCGCGCCGCACTGCTCCG'
seqGal1_TATA = 'ATATATAAA'
seqGal10_TATA = 'TTATATG'

seqGal10 = 'TCAGGAAAATCTGTAGACAATCTTGGACCCGTAAGTTTCACCGTTTTTCAAGGTTACACAATCTTTCCAGTTCTCTTGATTGATAGCATCAATGTATCTACCAGGCTCAATTGCAAAACCTTGTCTTGCTTCGTAACCAGCAGACAAGAAATCACCGGTATAAAATTGATAAGTTGGCTCTGTACTTAAAACTTCTAATGTAATATTGGAATCGGGATGAAAAGCCTTGACAATAAGCGTCAATTCATTGTTTAGAGTATTGATTTGACTTGGCTTAGCATTTTCATCCACCACAAAACAACAATCAAACTGGGGATTTTTGGGGCCTAAGACCGTTGGCTTTGTAGAGTTAAAGGTAGCAATTTCTCTATCGACGATATTACCCGTAGGAATCATGTTTTTGTCGACATCAACAGATTTTTTTGAACGCACCATAATCTCCGTACCCTCAATAGTGTCTCCATATGGCTTGTTCAGATTGAAATAACTATGATTTGTTAAATTTATTGGCGTCGCTTCACCAGCAGTCAATTTACCTTTATATACCATTTCCAAACTTTTTTGGGCAACGTTCACAGTATACTGTATGGTTACCAATAGATCACCTGGAAATTCGGTGTCCTTCTCATTATCTATCAGCATGTACTCGGCGGTAAAAACATCCTTTGAAGGATTTTGAATGATGGGTCCCAAAAATCTTTTTCTGTGGAAAGAACCGATACTACTATGATTCGCATTAACGCCGTTATTAACGGTTAACTGATAGTCTTTGTTGCATAAACTAAACTTACCCTTCGAAATACGATTAGCATACCTGCCGATCGTGGCGCCTATATAAGCACTATCAGGATTCAAATACCCTTCCTCATTTTCATAGCCAAGAACAACTGATTGTCCGTTCACTTTCAGGTCAACAATGCTGGCGCCCAAATTGGCAAACGTGGCTTGAAATCTGGTGCCGGCACCAATAGTCACAAATCTTGCGTCATAACGCATATCTTCAGCGGAAAATCTGGCCTCGACACCCCTTAACTGGTAACCAAAAGGATTCTCAGTAGTCCATTTCCATAAATCCTTGCAGGAGTCTTCAACCTGCAACTCGGTCTGCCATTTCAGTTCGCGTTTGGCCCTATCTGGTTTAGCCGTCAAGTTCAAAACATCACCTGCTCTTCTGCCCGTAACTTTGTATGGAAGATCAATACCAGAAGCTTTGCAGAATGCATGATAAACTTCAAAAACTGTAGAACCTTTACCGGAACCCAAGTTCCACTCACGACACAAACCTTCATTTTCATTGTAGGCCTCTAGGTATTGCAGGGCTGCAATATGACCTTTTGCTAGATCAACTACGTGGATATAATCCCTGATCGGGGTACCATCTCTGGAATCATAATCGTCTCCGAAGATGTAAAGCTTCTCGCGCCTACCAACAGCTACTTGAGCCATATATGGCAACAAATTGTTTGGTATACCTAGCGGATCTTCTCCGATTAATCCAGAGGGATGTGCGCCAATTGGGTTAAAATAACGCAAGATAGCAAACTTCCAACTTTTTTTGTCGCTATTGTAAAGATCATTCAAGATATTCTCAATGGCGTATTTCGTATGACCATACGGATTAGTAGGCCCTAAGGGACATTCTTCTGGGATAGGAATCATATTTGGGAATCTCGTAGCATCACCATAGACAGTAGCAGAAGATGAAAAAACAAATTTGGAAACGTTGTATTGTTGCATTAACTCTAATAAAACGACAGTTCCCAAAATGTTATTGTGATAGTATCTCAGCGGGATTTGTGTAGATTCACCTACAGCCTTTAAACCAGCAAAGTGAATTACCGAATCAATTTTATATTCTTTGAAAACCTTTTCCAGACCTTTTCGGTCACACAAATCAACCTCATAGAAGGGAATGTGATGCTTGGTCAAGACCTCTAACCTGGCTACAGAATCATAAGTTGAATTCGACAGGTTATCAGCAACAACACAGTCATATCCATTCTCAATTAGCTCTACCACAGTGTGTGAACCAATGTATCCAGCACCACCTGTAACCAAAACAATTTTAGAAGTACTTTCACTTTGTAACTGAGCTGTCAT'
seqGal1 = 'ATGACTAAATCTCATTCAGAAGAAGTGATTGTACCTGAGTTCAATTCTAGCGCAAAGGAATTACCAAGACCATTGGCCGAAAAGTGCCCGAGCATAATTAAGAAATTTATAAGCGCTTATGATGCTAAACCGGATTTTGTTGCTAGATCGCCTGGTAGAGTCAATCTAATTGGTGAACATATTGATTATTGTGACTTCTCGGTTTTACCTTTAGCTATTGATTTTGATATGCTTTGCGCCGTCAAAGTTTTGAACGAGAAAAATCCATCCATTACCTTAATAAATGCTGATCCCAAATTTGCTCAAAGGAAGTTCGATTTGCCGTTGGACGGTTCTTATGTCACAATTGATCCTTCTGTGTCGGACTGGTCTAATTACTTTAAATGTGGTCTCCATGTTGCTCACTCTTTTCTAAAGAAACTTGCACCGGAAAGGTTTGCCAGTGCTCCTCTGGCCGGGCTGCAAGTCTTCTGTGAGGGTGATGTACCAACTGGCAGTGGATTGTCTTCTTCGGCCGCATTCATTTGTGCCGTTGCTTTAGCTGTTGTTAAAGCGAATATGGGCCCTGGTTATCATATGTCCAAGCAAAATTTAATGCGTATTACGGTCGTTGCAGAACATTATGTTGGTGTTAACAATGGCGGTATGGATCAGGCTGCCTCTGTTTGCGGTGAGGAAGATCATGCTCTATACGTTGAGTTCAAACCGCAGTTGAAGGCTACTCCGTTTAAATTTCCGCAATTAAAAAACCATGAAATTAGCTTTGTTATTGCGAACACCCTTGTTGTATCTAACAAGTTTGAAACCGCCCCAACCAACTATAATTTAAGAGTGGTAGAAGTCACTACAGCTGCAAATGTTTTAGCTGCCACGTACGGTGTTGTTTTACTTTCTGGAAAAGAAGGATCGAGCACGAATAAAGGTAATCTAAGAGATTTCATGAACGTTTATTATGCCAGATATCACAACATTTCCACACCCTGGAACGGCGATATTGAATCCGGCATCGAACGGTTAACAAAGATGCTAGTACTAGTTGAAGAGTCTCTCGCCAATAAGAAACAGGGCTTTAGTGTTGACGATGTCGCACAATCCTTGAATTGTTCTCGCGAAGAATTCACAAGAGACTACTTAACAACATCTCCAGTGAGATTTCAAGTCTTAAAGCTATATCAGAGGGCTAAGCATGTGTATTCTGAATCTTTAAGAGTCTTGAAGGCTGTGAAATTAATGACTACAGCGAGCTTTACTGCCGACGAAGACTTTTTCAAGCAATTTGGTGCCTTGATGAACGAGTCTCAAGCTTCTTGCGATAAACTTTACGAATGTTCTTGTCCAGAGATTGACAAAATTTGTTCCATTGCTTTGTCAAATGGATCATATGGTTCCCGTTTGACCGGAGCTGGCTGGGGTGGTTGTACTGTTCACTTGGTTCCAGGGGGCCCAAATGGCAACATAGAAAAGGTAAAAGAAGCCCTTGCCAATGAGTTCTACAAGGTCAAGTACCCTAAGATCACTGATGCTGAGCTAGAAAATGCTATCATCGTCTCTAAACCAGCATTGGGCAGCTGTCTATATGAATTATAA'

seqPromGal10GAL1 = 'TCGGTCTGCCATTTCAGTTCGCGTTTGGCCCTATCTGGTTTAGCCGTCAAGTTCAAAACATCACCTGCTCTTCTGCCCGTAACTTTGTATGGAAGATCAATACCAGAAGCTTTGCAGAATGCATGATAAACTTCAAAAACTGTAGAACCTTTACCGGAACCCAAGTTCCACTCACGACACAAACCTTCATTTTCATTGTAGGCCTCTAGGTATTGCAGGGCTGCAATATGACCTTTTGCTAGATCAACTACGTGGATATAATCCCTGATCGGGGTACCATCTCTGGAATCATAATCGTCTCCGAAGATGTAAAGCTTCTCGCGCCTACCAACAGCTACTTGAGCCATATATGGCAACAAATTGTTTGGTATACCTAGCGGATCTTCTCCGATTAATCCAGAGGGATGTGCGCCAATTGGGTTAAAATAACGCAAGATAGCAAACTTCCAACTTTTTTTGTCGCTATTGTAAAGATCATTCAAGATATTCTCAATGGCGTATTTCGTATGACCATACGGATTAGTAGGCCCTAAGGGACATTCTTCTGGGATAGGAATCATATTTGGGAATCTCGTAGCATCACCATAGACAGTAGCAGAAGATGAAAAAACAAATTTGGAAACGTTGTATTGTTGCATTAACTCTAATAAAACGACAGTTCCCAAAATGTTATTGTGATAGTATCTCAGCGGGATTTGTGTAGATTCACCTACAGCCTTTAAACCAGCAAAGTGAATTACCGAATCAATTTTATATTCTTTGAAAACCTTTTCCAGACCTTTTCGGTCACACAAATCAACCTCATAGAAGGGAATGTGATGCTTGGTCAAGACCTCTAACCTGGCTACAGAATCATAAGTTGAATTCGACAGGTTATCAGCAACAACACAGTCATATCCATTCTCAATTAGCTCTACCACAGTGTGTGAACCAATGTATCCAGCACCACCTGTAACCAAAACAATTTTAGAAGTACTTTCACTTTGTAACTGAGCTGTCATTTATATTGAATTTTCAAAAATTCTTACTTTTTTTTTGGATGGACGCAAAGAAGTTTAATAATCATATTACATGGCATTACCACCATATACATATCCATATCTAATCTTACTTATATGTTGTGGAAATGTAAAGAGCCCCATTATCTTAGCCTAAAAAAACCTTCTCTTTGGAACTTTCAGTAATACGCTTAACTGCTCATTGCTATATTGAAGTACGGATTAGAAGCCGCCGAGCGGGCGACAGCCCTCCGACGGAAGACTCTCCTCCGTGCGTCCTCGTCTTCACCGGTCGCGTTCCTGAAACGCAGATGTGCCTCGCGCCGCACTGCTCCGAACAATAAAGATTCTACAATACTAGCTTTTATGGTTATGAAGAGGAAAAATTGGCAGTAACCTGGCCCCACAAACCTTCAAATTAACGAATCAAATTAACAACCATAGGATGATAATGCGATTAGTTTTTTAGCCTTATTTCTGGGGTAATTAATCAGCGAAGCGATGATTTTTGATCTATTAACAGATATATAAATGGAAAAGCTGCATAACCACTTTAACTAATACTTTCAACATTTTCAGTTTGTATTACTTCTTATTCAAATGTCATAAAAGTATCAACAAAAAATTGTTAATATACCTCTATACTTTAACGTCAAGGAGAAAAAACTATAATGACTAAATCTCATTCAGAAGAAGTGATTGTACCTGAGTTCAATTCTAGCGCAAAGGAATTACCAAGACCATTGGCCGAAAAGTGCCCGAGCATAATTAAGAAATTTATAAGCGCTTATGATGCTAAACCGGATTTTGTTGCTAGATCGCCTGGTAGAGTCAATCTAATTGGTGAACATATTGATTATTGTGACTTCTCGGTTTTACCTTTAGCTATTGATTTTGATATGCTTTGCGCCGTCAAAGTTTTGAACGAGAAAAATCCATCCATTACCTTAATAAATGCTGATCCCAAATTTGCTCAAAGGAAGTTCGATTTGCCGTTGGACGG'


# In[11]:


## function to get the LLR data out of the data of first reconstitution 
def read_stats(filename, seq, length=0, range=None, motifs=None, motifs2=None, askstrand=''):
    if 'barcode01' in filename:
        if askstrand == '+':
            strands = [0]
        elif askstrand == '-':
            strands = [1]
        else:
            strands = [0,1]
    else:
        if askstrand == '+':
            strands = [1]
        elif askstrand == '-':
            strands = [0]
        else:
            strands = [1, 0]

    result = []
    strand = []
    for block in strands:
        df = pd.read_hdf(filename, key=f'Statistic_Blocks/Block_{block}/block_stats')
        for id in tqdm.tqdm(df['read_id'].unique(), desc=f'Reading strand: {block}'):
            mod = np.zeros((len(seq)))
            df_id = df[df['read_id'] == id]
            mod[np.asarray(df_id['pos'])] = np.asarray(df_id['stat'])
            if range is None:
                if (df_id['pos'].max() - df_id['pos'].min()) > length and df_id['pos'].min():
                    result.append(mod)
                    strand.append(block)
            else:
                if (df_id['pos'].max() >= range[1] and df_id['pos'].min()) <= range[0]:
                    result.append(mod)
                    strand.append(block)
    result = np.asarray(result)

    if motifs is not None:
        filter = np.zeros(len(seq))
        for motif in motifs:
            index = find_motif(seq, motif, format='index')
            offset = [i for i, c in enumerate(motif) if c.isupper()]
            filter[index + offset[0]] = 1
        result = np.multiply(result, filter)
    
    if motifs2 is not None:
        filter2 = np.ones(len(seq))
        for motif2 in motifs2:
            index = find_motif(seq, motif2, format='index')
            offset = [i for i, c in enumerate(motif2) if c.isupper()]
            if len(index + offset[0]) != 0:
                filter2[index + offset[0]] = 0
        result = np.multiply(result, filter2) 
            
    return result, np.asarray(strand)



if __name__ == '__main__':
    barcode=1
    barcode1 = 1
    barcode2 = 2
    barcode3 = 3

    directory = r'/DATA/lenstra_lab/c.blom/Tests/Test_sequence_aligning_reconstitution1'
    #directory = r'/DATA/lenstra_lab/c.blom/Tests/Reconsitution_2ndExp_methylation'

    experiment = directory+r'/sample_descriptionTestbarcode3.xlsx'
    #experiment = directory+r'/sample_descriptionTestbarcodes1_4.xlsx'

    # tombo_filename = fr'/media/noort/Data/users/noort/20220816_barcode{barcode:02d}_selection/read_stats.MOD.tombo.per_read_stats'
    tombo_filename1 = directory + fr'/barcode{barcode1:02d}.read_stats.MOD.tombo.per_read_stats'
    tombo_filename2 = directory + fr'/barcode{barcode2:02d}.read_stats.MOD.tombo.per_read_stats'
    tombo_filename3 = directory + fr'/barcode{barcode3:02d}.read_stats.MOD.tombo.per_read_stats'
    fasta = directory+rf'/LinearGalInplasmid.fasta'

    mods = ['5mC', '6mA']
    # fast5_filename = directory+r'/fast5_pass/barcode{barcode:02d}/tombo/0/0b1f2bf3-2287-4af9-8dd4-5c26d11269fb'

    df = pd.read_excel(experiment, index_col='Barcode')
    label = df.at[barcode, 'Enzymes']
    print(df)


    for mod in mods[:-1]:
        files=[tombo_filename1.replace('MOD', mod)] #, tombo_filename2.replace('MOD', mod), tombo_filename4.replace('MOD', mod)]
        files=np.array(files, dtype='object')
        barcodes=['barcode01','barcode02', 'barcode03']

        askstrand='+'

        llrRec1Barcode1plus, new_llrRec1Barcode1plus, strandsRec1Barcode1plus = read_files(tombo_filename1.replace('MOD', mod), fasta, askstrand)
        llrRec1Barcode2plus, new_llrRec1Barcode2plus, strandsRec1Barcode2plus = read_files(tombo_filename2.replace('MOD', mod), fasta, askstrand)
        llrRec1Barcode3plus, new_llrRec1Barcode3plus, strandsRec1Barcode3plus = read_files(tombo_filename3.replace('MOD', mod), fasta, askstrand)
        
        askstrand='-'

        llrRec1Barcode1min, new_llrRec1Barcode1min, strandsRec1Barcode1min = read_files(tombo_filename1.replace('MOD', mod), fasta, askstrand)
        llrRec1Barcode2min, new_llrRec1Barcode2min, strandsRec1Barcode2min = read_files(tombo_filename2.replace('MOD', mod), fasta, askstrand)
        llrRec1Barcode3min, new_llrRec1Barcode3min, strandsRec1Barcode3min = read_files(tombo_filename3.replace('MOD', mod), fasta, askstrand)


# In[12]:


## function to get the LLR data out of the data of second reconstitution 
def read_stats(filename, seq, length=0, range=None, motifs=None, motifs2=None, askstrand=''):
    if askstrand == '+':
        strands = [1]
    elif askstrand == '-':
        strands = [0]
    else:
        strands = [1, 0]

    result = []
    strand = []
    for block in strands:
        df = pd.read_hdf(filename, key=f'Statistic_Blocks/Block_{block}/block_stats')
        for id in tqdm.tqdm(df['read_id'].unique(), desc=f'Reading strand: {block}'):
            mod = np.zeros((len(seq)))
            df_id = df[df['read_id'] == id]
            mod[np.asarray(df_id['pos'])] = np.asarray(df_id['stat'])
            if range is None:
                if (df_id['pos'].max() - df_id['pos'].min()) > length and df_id['pos'].min():
                    result.append(mod)
                    strand.append(block)
            else:
                if (df_id['pos'].max() >= range[1] and df_id['pos'].min()) <= range[0]:
                    result.append(mod)
                    strand.append(block)
    result = np.asarray(result)

    if motifs is not None:
        filter = np.zeros(len(seq))
        for motif in motifs:
            index = find_motif(seq, motif, format='index')
            offset = [i for i, c in enumerate(motif) if c.isupper()]
            filter[index + offset[0]] = 1
        result = np.multiply(result, filter)
        
    if motifs2 is not None:
        filter2 = np.ones(len(seq))
        for motif2 in motifs2:
            index = find_motif(seq, motif2, format='index')
            offset = [i for i, c in enumerate(motif2) if c.isupper()]
            if len(index + offset[0]) != 0:
                filter2[index + offset[0]] = 0
        result = np.multiply(result, filter2) 
        
    return result, np.asarray(strand)

if __name__ == '__main__':
    barcode=1
    barcode1 = 1
    barcode2 = 2
    barcode3 = 3

    #directory = r'/DATA/lenstra_lab/c.blom/Tests/Test_sequence_aligning_reconstitution1'
    directory = r'/DATA/lenstra_lab/c.blom/Tests/Reconsitution_2ndExp_methylation'


    experiment = directory+r'/sample_descriptionTestbarcode3.xlsx'
    #experiment = directory+r'/sample_descriptionTestbarcodes1_4.xlsx'

    # tombo_filename = fr'/media/noort/Data/users/noort/20220816_barcode{barcode:02d}_selection/read_stats.MOD.tombo.per_read_stats'
    tombo_filename1 = directory + fr'/barcode{barcode1:02d}.read_stats.MOD.tombo.per_read_stats'
    tombo_filename2 = directory + fr'/barcode{barcode2:02d}.read_stats.MOD.tombo.per_read_stats'
    tombo_filename3 = directory + fr'/barcode{barcode3:02d}.read_stats.MOD.tombo.per_read_stats'
    fasta = directory+rf'/LinearGalInplasmid.fasta'

    mods = ['5mC', '6mA']
    # fast5_filename = directory+r'/fast5_pass/barcode{barcode:02d}/tombo/0/0b1f2bf3-2287-4af9-8dd4-5c26d11269fb'

    df = pd.read_excel(experiment, index_col='Barcode')
    label = df.at[barcode, 'Enzymes']
    print(df)

    for mod in mods[:-1]:
        files=[tombo_filename1.replace('MOD', mod)] #, tombo_filename2.replace('MOD', mod), tombo_filename4.replace('MOD', mod)]
        files=np.array(files, dtype='object')
        barcodes=['barcode01','barcode02', 'barcode03']

        askstrand='+'

        llrRec2Barcode1plus, new_llrRec2Barcode1plus, strandsRec2Barcode1plus = read_files(tombo_filename1.replace('MOD', mod), fasta, askstrand)
        llrRec2Barcode2plus, new_llrRec2Barcode2plus, strandsRec2Barcode2plus = read_files(tombo_filename2.replace('MOD', mod), fasta, askstrand)
        llrRec2Barcode3plus, new_llrRec2Barcode3plus, strandsRec2Barcode3plus = read_files(tombo_filename3.replace('MOD', mod), fasta, askstrand)
        
        askstrand='-'

        llrRec2Barcode1min, new_llrRec2Barcode1min, strandsRec2Barcode1min = read_files(tombo_filename1.replace('MOD', mod), fasta, askstrand)
        llrRec2Barcode2min, new_llrRec2Barcode2min, strandsRec2Barcode2min = read_files(tombo_filename2.replace('MOD', mod), fasta, askstrand)
        llrRec2Barcode3min, new_llrRec2Barcode3min, strandsRec2Barcode3min = read_files(tombo_filename3.replace('MOD', mod), fasta, askstrand)


# In[19]:


## Testing how many principical components you must have to take into account most data points
LLRfilename = llrRec1Barcode1plus
motif = seqProm

seq = read_fasta(fasta).upper()

index = find_motif(seq, motif, format='index')[0]
seq = find_motif(seq, motif, format='string')
xlabels = [s + f'\n:\n{i}' if i % 10 == 0 else s for i, s in enumerate(seq)]

plotrange= [index, index+len(motif)]
subllr = LLRfilename[:, plotrange[0]:plotrange[1]]
    
scaler = MaxAbsScaler()
X_scaled = scaler.fit_transform(subllr)
    
pca= PCA(n_components = 2, random_state = 42)
X_principal = pca.fit_transform(X_scaled)
a=pca.explained_variance_ratio_.cumsum()
print(a)


# In[22]:


new_indexplus, new_indexmin = find_index_methylation(tombo_filename1.replace('MOD', mod), fasta)
outputdir = r'/DATA/lenstra_lab/c.blom/Tests/Figures/'


# In[108]:


# figure 3c: GAL locus without WCE, + and - strand, no clustering, 1st reconstitution
llr = np.concatenate((llrRec1Barcode3min, llrRec1Barcode3plus))
strands = np.append(strandsRec1Barcode3min, strandsRec1Barcode3plus)

#plot_llr(llrRec1Barcode3plus, strandsRec1Barcode3plus, fasta, seqGAL, new_indexplus, 'sum')
plot_llr(llr, strands, fasta, seqGAL, new_indexplus, outputdir,'figure3c.pdf', 'sum')


# In[64]:


n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(llrRec1Barcode3plus, fasta, seqGAL, new_indexplus, outputdir,'figure3dSup1.pdf')
dictClusterNoWCE = split_clusters(llrRec1Barcode3plus, n_clusterfinal, predicted_labels)
dictClusterNoWCEnotEmpty = dictClusterNoWCE['Cluster0']
# figure 3d supplementary: GAL locus without WCE, with empty reads, K-means over whole GAL locus, 1st reconstitution
plot_llr(llrRec1Barcode3plus, '', fasta, seqGAL, new_indexplus,  outputdir,'figure3dSup2.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 3d: reclustering GAL locus without WCE, without empty reads, K-means over whole GAL locus, 1st reconstitution
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterNoWCE['Cluster0'], fasta, seqGAL, new_indexplus,  outputdir,'figure3d1.pdf')

plot_llr(dictClusterNoWCE['Cluster0'], '', fasta, seqGAL, new_indexplus,  outputdir,'figure3d2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2)


# In[65]:


# figure 3e supplementary: Hierarchical clustering on GAL locus with empty reads of No WCE reconstitution
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_agglomerative_determinant(llrRec1Barcode3plus, fasta, seqGAL, new_indexplus, outputdir,'figure3eSup1.pdf')
new_direction = plot_dendrogram(llrRec1Barcode3plus, pcadf, outputdir,'figure3eSup2.pdf')
plot_llr(new_direction, '', fasta, seqGAL, new_indexplus, outputdir,'figure3eSup3.pdf', method = 'dendrogram', predicted_labels = None, ncluster = None)

# figure 3e supplementary: Hierarchical clustering on GAL locus without empty reads of No WCE reconstitution
dictClusterNoWCEHierarchical = split_clusters(llrRec1Barcode3plus, n_clusterfinal, predicted_labels)
dictClusterNoWCEnotEmptyHierarchical = dictClusterNoWCEHierarchical['Cluster0']

n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_agglomerative_determinant(dictClusterNoWCEHierarchical['Cluster0'], fasta, seqGAL, new_indexplus, outputdir,'figure3e1.pdf')
new_direction = plot_dendrogram(dictClusterNoWCEHierarchical['Cluster0'], pcadf, outputdir,'figure3e2.pdf')
plot_llr(new_direction, '', fasta, seqGAL, new_indexplus, outputdir,'figure3e3.pdf', 'dendrogram', predicted_labels = None, ncluster = None)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100
cluster3 = np.count_nonzero(predicted_labels == 3)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2, cluster3)


# In[70]:


n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(llrRec1Barcode1plus, fasta, seqGAL, new_indexplus, outputdir,'figure4aSup1.pdf')
dictClusterGALWCE = split_clusters(llrRec1Barcode1plus, n_clusterfinal, predicted_labels)
dictClusterGALWCEnotEmpty = dictClusterGALWCE['Cluster1']
# figure 4a supplementary: GAL locus with GAL WCE, with empty reads, K-means over whole GAL locus, 1st reconstitution
plot_llr(llrRec1Barcode1plus, '', fasta, seqGAL, new_indexplus, outputdir,'figure4aSup2.pdf','clusters', predicted_labels, n_clusterfinal)

# figure 4a: reclustering GAL locus with GAL, without empty reads, K-means over whole GAL locus, 1st reconstitution
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterGALWCE['Cluster1'], fasta, seqGAL, new_indexplus, outputdir,'figure4a1.pdf')

plot_llr(dictClusterGALWCE['Cluster1'], '', fasta, seqGAL, new_indexplus, outputdir,'figure4a2.pdf','clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2)


# In[71]:


# figure 4b supplementary: Hierarchical clustering on GAL locus with empty reads of GAL WCE reconstitution
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_agglomerative_determinant(llrRec1Barcode1plus, fasta, seqGAL, new_indexplus, outputdir,'figure4bSup1.pdf')
new_direction = plot_dendrogram(llrRec1Barcode1plus, pcadf, outputdir,'figure4bSup2.pdf')
plot_llr(new_direction, '', fasta, seqGAL, new_indexplus, outputdir,'figure4bSup3.pdf', method = 'dendrogram', predicted_labels = None, ncluster = None)

# figure 4b supplementary: Hierarchical clustering on GAL locus without empty reads of GAL WCE reconstitution
dictClusterGALWCEHierarchical = split_clusters(llrRec1Barcode1plus, n_clusterfinal, predicted_labels)
dictClusterGALWCEnotEmptyHierarchical = dictClusterGALWCEHierarchical['Cluster0']

n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_agglomerative_determinant(dictClusterGALWCEHierarchical['Cluster0'], fasta, seqGAL, new_indexplus, outputdir,'figure4b1.pdf')
new_direction = plot_dendrogram(dictClusterGALWCEHierarchical['Cluster0'], pcadf, outputdir,'figure4b2.pdf')
plot_llr(new_direction, '', fasta, seqGAL, new_indexplus, outputdir,'figure4b3.pdf', 'dendrogram', predicted_labels = None, ncluster = None)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1)


# In[74]:


n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(llrRec1Barcode2plus, fasta, seqGAL, new_indexplus, outputdir,'figure4dSup1.pdf')
dictClusterRAFWCE = split_clusters(llrRec1Barcode2plus, n_clusterfinal, predicted_labels)
dictClusterRAFWCEnotEmpty = dictClusterRAFWCE['Cluster1']
# figure 4d supplementary: GAL locus with GAL WCE, with empty reads, K-means over whole GAL locus, 1st reconstitution
plot_llr(llrRec1Barcode2plus, '', fasta, seqGAL, new_indexplus, outputdir,'figure4dSup2.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 4d: reclustering GAL locus with GAL, without empty reads, K-means over whole GAL locus, 1st reconstitution
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterRAFWCE['Cluster1'], fasta, seqGAL, new_indexplus, outputdir,'figure4d1.pdf')

plot_llr(dictClusterRAFWCE['Cluster1'], '', fasta, seqGAL, new_indexplus, outputdir,'figure4d2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2)


# In[75]:


# figure 4e supplementary: Hierarchical clustering on GAL locus with empty reads of RAF WCE reconstitution
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_agglomerative_determinant(llrRec1Barcode2plus, fasta, seqGAL, new_indexplus, outputdir,'figure4eSup1.pdf')
new_direction = plot_dendrogram(llrRec1Barcode2plus, pcadf, outputdir,'figure4eSup2.pdf')
plot_llr(new_direction, '', fasta, seqGAL, new_indexplus, outputdir,'figure4eSup3.pdf', method = 'dendrogram', predicted_labels = None, ncluster = None)

# figure 4e: Hierarchical clustering on GAL locus without empty reads of RAF WCE reconstitution
dictClusterRAFWCEHierarchical = split_clusters(llrRec1Barcode2plus, n_clusterfinal, predicted_labels)
dictClusterRAFWCEnotEmptyHierarchical = dictClusterRAFWCEHierarchical['Cluster1']

n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_agglomerative_determinant(dictClusterRAFWCEHierarchical['Cluster1'], fasta, seqGAL, new_indexplus, outputdir,'figure4e1.pdf')
new_direction = plot_dendrogram(dictClusterRAFWCEHierarchical['Cluster1'], pcadf, outputdir,'figure4e2.pdf')
plot_llr(new_direction, '', fasta, seqGAL, new_indexplus, outputdir,'figure4e3.pdf', 'dendrogram', predicted_labels = None, ncluster = None)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1)


# In[130]:


# block 21
a = dictClusterNoWCEnotEmpty
b = dictClusterGALWCEnotEmpty
c = dictClusterRAFWCEnotEmpty

motif = seqGAL
plotrange, withinrange = llr_in_range(new_indexplus, fasta, motif)

DifGal = SumPerSiteNewVersion(b)-SumPerSiteNewVersion(a)
DifRaf = SumPerSiteNewVersion(c)-SumPerSiteNewVersion(a)

DifGalwithRAF = SumPerSiteNewVersion(b)-SumPerSiteNewVersion(c)

# figure 4c: difference GAL WCE and no WCE
namefig = 'figure4c.pdf'
plt.figure(figsize = (14,4) ,dpi=200)
markerline, stemlines, baseline = plt.stem(DifGal, use_line_collection = True)
plt.setp(stemlines, 'linewidth', 2)

plt.xlim(plotrange[0],plotrange[1])
plt.xlabel('i')
plt.ylim(-1,1)
plt.ylabel('\u0394 average LLR (GAL - No WCE) ')
plt.savefig(outputdir + namefig, dpi=1200, transparent=True)
plt.show()

# figure 4f: difference RAF WCE and no WCE
namefig = 'figure4f.pdf'
plt.figure(figsize = (14,4) ,dpi=200)
markerline, stemlines, baseline = plt.stem(DifRaf, use_line_collection = True)
plt.setp(stemlines, 'linewidth', 2)

plt.ylim(-1,1)
plt.ylabel('\u0394 average LLR (RAF - No WCE)')
plt.xlim(plotrange[0],plotrange[1])
plt.xlabel('i')
plt.savefig(outputdir + namefig, dpi=1200, transparent=True)
plt.show()

# figure 
plt.figure(figsize = (14,4) ,dpi=200)
markerline, stemlines, baseline = plt.stem(DifGalwithRAF, use_line_collection = True)
plt.setp(stemlines, 'linewidth', 2)

plt.xlim(plotrange[0],plotrange[1])
plt.xlabel('i')
plt.ylim(-1,1)
plt.ylabel('\u0394 average LLR (GAL - RAF) ')
plt.show()


# In[119]:


# figure 5a: reclustering on Promoter of No WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterNoWCEnotEmpty, fasta, seqProm, new_indexplus, outputdir,'figure5a1.pdf')
plot_llr(dictClusterNoWCEnotEmpty, '', fasta, seqProm, new_indexplus, outputdir, 'figure5a2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2)

# supplementary figure 5a:
plot_llr(dictClusterNoWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5aSup.pdf','clusters', predicted_labels, n_clusterfinal)

# figure 5b: reclustering on Promoter of GAL WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterGALWCEnotEmpty, fasta, seqProm, new_indexplus, outputdir,'figure5b1.pdf')
plot_llr(dictClusterGALWCEnotEmpty, '', fasta, seqProm, new_indexplus, outputdir,'figure5b2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1)

# supplementary figure 5b:
plot_llr(dictClusterGALWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5bSup.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 5c: reclustering on Promoter of RAF WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterRAFWCEnotEmpty, fasta, seqProm, new_indexplus, outputdir,'figure5c1.pdf')
plot_llr(dictClusterRAFWCEnotEmpty, '', fasta, seqProm, new_indexplus, outputdir,'figure5c2.pdf','clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1)

# supplementary figure 5c:
plot_llr(dictClusterRAFWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5cSup.pdf', 'clusters', predicted_labels, n_clusterfinal)


# figure 5d: Sumplots to compare different conditions
plotSumPerSiteComparingConditions(dictClusterGALWCEnotEmpty, dictClusterRAFWCEnotEmpty, dictClusterNoWCEnotEmpty, fasta, new_indexplus, seqProm, outputdir,'figure5d.pdf')


# In[120]:


# figure 5e: reclustering on Promoter + Gal10 of No WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterNoWCEnotEmpty, fasta, seqProm_Gal10, new_indexplus, outputdir,'figure5e1.pdf')
plot_llr(dictClusterNoWCEnotEmpty, '', fasta, seqProm_Gal10, new_indexplus, outputdir,'figure5e2.pdf','clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2)

# supplementary figure 5e:
plot_llr(dictClusterNoWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5eSup.pdf','clusters', predicted_labels, n_clusterfinal)

# figure 5f: reclustering on Promoter + Gal10 of GAL WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterGALWCEnotEmpty, fasta, seqProm_Gal10, new_indexplus, outputdir,'figure5f1.pdf')
plot_llr(dictClusterGALWCEnotEmpty, '', fasta, seqProm_Gal10, new_indexplus, outputdir,'figure5f2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1)

# supplementary figure 5f:
plot_llr(dictClusterGALWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5fSup.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 5g: reclustering on Promoter + Gal10 of RAF WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterRAFWCEnotEmpty, fasta, seqProm_Gal10, new_indexplus, outputdir,'figure5g1.pdf')
plot_llr(dictClusterRAFWCEnotEmpty, '', fasta, seqProm_Gal10, new_indexplus, outputdir,'figure5g2.pdf','clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2)

# supplementary figure 5g:
plot_llr(dictClusterRAFWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5gSup.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 5h: Sumplots to compare different conditions
plotSumPerSiteComparingConditions(dictClusterGALWCEnotEmpty, dictClusterRAFWCEnotEmpty, dictClusterNoWCEnotEmpty, fasta, new_indexplus, seqProm_Gal10, outputdir,'figure5h.pdf',)


# In[121]:


# figure 5i: reclustering on Promoter + Gal1 of No WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterNoWCEnotEmpty, fasta, seqProm_Gal1, new_indexplus, outputdir,'figure5i1.pdf')
plot_llr(dictClusterNoWCEnotEmpty, '', fasta, seqProm_Gal1, new_indexplus, outputdir,'figure5i2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100
cluster3 = np.count_nonzero(predicted_labels == 3)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2, cluster3)

# supplementary figure 5i:
plot_llr(dictClusterNoWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5iSup.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 5j: reclustering on Promoter + Gal1 of GAL WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterGALWCEnotEmpty, fasta, seqProm_Gal1, new_indexplus, outputdir,'figure5j1.pdf')
plot_llr(dictClusterGALWCEnotEmpty, '', fasta, seqProm_Gal1, new_indexplus, outputdir,'figure5j2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2)

# supplementary figure 5j:
plot_llr(dictClusterGALWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5jSup.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 5k: reclustering on Promoter + Gal1 of RAF WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterRAFWCEnotEmpty, fasta, seqProm_Gal1, new_indexplus, outputdir,'figure5k1.pdf')
plot_llr(dictClusterRAFWCEnotEmpty, '', fasta, seqProm_Gal1, new_indexplus, outputdir,'figure5k2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2)

# supplementary figure 5k:
plot_llr(dictClusterRAFWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5kSup.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 5l: Sumplots to compare different conditions
plotSumPerSiteComparingConditions(dictClusterGALWCEnotEmpty, dictClusterRAFWCEnotEmpty, dictClusterNoWCEnotEmpty, fasta, new_indexplus, seqProm_Gal1, outputdir,'figure5l.pdf')


# In[138]:


# figure 5m: reclustering on GAL10 of No WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterNoWCEnotEmpty, fasta, seqGal10, new_indexplus, outputdir,'figure5m1.pdf')
plot_llr(dictClusterNoWCEnotEmpty, '', fasta, seqGal10, new_indexplus, outputdir,'figure5m2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2)

NoWCEGAL10 = split_clusters(dictClusterNoWCEnotEmpty, n_clusterfinal, predicted_labels)

NoWCEGAL10num0 = NoWCEGAL10['Cluster0']
NoWCEGAL10num1 = NoWCEGAL10['Cluster1']
NoWCEGAL10num2 = NoWCEGAL10['Cluster2']

# supplementary figure 5m:
plot_llr(dictClusterNoWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5mSup.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 5n: reclustering on GAL10 of GAL WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterGALWCEnotEmpty, fasta, seqGal10, new_indexplus, outputdir,'figure5n1.pdf')
plot_llr(dictClusterGALWCEnotEmpty, '', fasta, seqGal10, new_indexplus, outputdir,'figure5n2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2)

GALWCEGAL10 = split_clusters(dictClusterGALWCEnotEmpty, n_clusterfinal, predicted_labels)

GALWCEGAL10num0 = GALWCEGAL10['Cluster0']
GALWCEGAL10num1 = GALWCEGAL10['Cluster1']
GALWCEGAL10num2 = GALWCEGAL10['Cluster2']

# supplementary figure 5n:
plot_llr(dictClusterGALWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5nSup.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 5o: reclustering on GAL10 of RAF WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterRAFWCEnotEmpty, fasta, seqGal10, new_indexplus, outputdir,'figure5o1.pdf')
plot_llr(dictClusterRAFWCEnotEmpty, '', fasta, seqGal10, new_indexplus, outputdir,'figure5o2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2)

RAFWCEGAL10 = split_clusters(dictClusterRAFWCEnotEmpty, n_clusterfinal, predicted_labels)

RAFWCEGAL10num0 = RAFWCEGAL10['Cluster0']
RAFWCEGAL10num1 = RAFWCEGAL10['Cluster1']
RAFWCEGAL10num2 = RAFWCEGAL10['Cluster2']

# supplementary figure 5o:
plot_llr(dictClusterRAFWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5oSup.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 5p: Sumplots to compare different conditions
plotSumPerSiteComparingConditions(dictClusterGALWCEnotEmpty, dictClusterRAFWCEnotEmpty, dictClusterNoWCEnotEmpty, fasta, new_indexplus, seqGal10, outputdir,'figure5p.pdf')


# In[140]:


# figure 5q: reclustering on GAL1 of No WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterNoWCEnotEmpty, fasta, seqGal1, new_indexplus, outputdir,'figure5q1.pdf')
plot_llr(dictClusterNoWCEnotEmpty, '', fasta, seqGal1, new_indexplus, outputdir,'figure5q2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100
cluster3 = np.count_nonzero(predicted_labels == 3)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2, cluster3)

NoWCEGAL1 = split_clusters(dictClusterNoWCEnotEmpty, n_clusterfinal, predicted_labels)

NoWCEGAL1num0 = NoWCEGAL1['Cluster0']
NoWCEGAL1num1 = NoWCEGAL1['Cluster1']
NoWCEGAL1num2 = NoWCEGAL1['Cluster2']
NoWCEGAL1num3 = NoWCEGAL1['Cluster3']

# supplementary figure 5q:
plot_llr(dictClusterNoWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5qSup.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 5r: reclustering on GAL1 of GAL WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterGALWCEnotEmpty, fasta, seqGal1, new_indexplus, outputdir,'figure5r1.pdf')
plot_llr(dictClusterGALWCEnotEmpty, '', fasta, seqGal1, new_indexplus, outputdir,'figure5r2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2)

GALWCEGAL1 = split_clusters(dictClusterGALWCEnotEmpty, n_clusterfinal, predicted_labels)

GALWCEGAL1num0 = GALWCEGAL1['Cluster0']
GALWCEGAL1num1 = GALWCEGAL1['Cluster1']
GALWCEGAL1num2 = GALWCEGAL1['Cluster2']

# supplementary figure 5r:
plot_llr(dictClusterGALWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5rSup.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 5s: reclustering on GAL1 of RAF WCE reconstitution using k-means
n_clusterfinal, silhouette_scores, pcadf, predicted_labels, plotrange = clustering_kmeans_determinantNewVersion(dictClusterRAFWCEnotEmpty, fasta, seqGal1, new_indexplus, outputdir,'figure5s1.pdf')
plot_llr(dictClusterRAFWCEnotEmpty, '', fasta, seqGal1, new_indexplus, outputdir,'figure5s2.pdf', 'clusters', predicted_labels, n_clusterfinal)

numberofdatapoints = np.size(predicted_labels)
cluster0 = np.count_nonzero(predicted_labels == 0)/numberofdatapoints*100
cluster1 = np.count_nonzero(predicted_labels == 1)/numberofdatapoints*100
cluster2 = np.count_nonzero(predicted_labels == 2)/numberofdatapoints*100
cluster3 = np.count_nonzero(predicted_labels == 3)/numberofdatapoints*100

print(numberofdatapoints, cluster0, cluster1, cluster2, cluster3)

RAFWCEGAL1 = split_clusters(dictClusterRAFWCEnotEmpty, n_clusterfinal, predicted_labels)

RAFWCEGAL1num0 = RAFWCEGAL1['Cluster0']
RAFWCEGAL1num1 = RAFWCEGAL1['Cluster1']
RAFWCEGAL1num2 = RAFWCEGAL1['Cluster2']
RAFWCEGAL1num3 = RAFWCEGAL1['Cluster3']

# supplementary figure 5s:
plot_llr(dictClusterRAFWCEnotEmpty, '', fasta, seqGAL, new_indexplus, outputdir,'figure5sSup.pdf', 'clusters', predicted_labels, n_clusterfinal)

# figure 5t: Sumplots to compare different conditions
plotSumPerSiteComparingConditions(dictClusterGALWCEnotEmpty, dictClusterRAFWCEnotEmpty, dictClusterNoWCEnotEmpty, fasta, new_indexplus, seqGal1, outputdir,'figure5t.pdf')


# In[199]:


# block 27
# to compare average with MNase-seq data from Ineke 
b = dictClusterGALWCEnotEmpty
c = dictClusterRAFWCEnotEmpty

motif = seqPromGal10GAL1
plotrange, withinrange = llr_in_range(new_indexplus, fasta, motif)

GALaverage = SumPerSiteNewVersion(b[:,withinrange])
RAFaverage = SumPerSiteNewVersion(c[:,withinrange])

# figure 
plt.figure(figsize = (14,4) ,dpi=200)
plt.step(withinrange, GALaverage, label = 'GAL WCE', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.step(withinrange, RAFaverage, label = 'RAF WCE', linewidth = 1, linestyle = '--') # scatter, s = 50, marker = 'o')
plt.ylabel('average LLR')
plt.ylim((-1,1))
plt.xlim(plotrange[0],plotrange[1])
plt.xlabel('i')
plt.legend(loc='upper right')
plt.grid(True)
plt.gca().invert_xaxis()

plt.savefig(outputdir + 'ComparewithMNase.pdf', dpi=1200, transparent=True)
plt.show()


# In[207]:


motif = seqGal10
plotrange, withinrange = llr_in_range(new_indexplus, fasta, motif)

# figure No WCE average clusters for GAL10
NoWCEGAL10num0average = SumPerSiteNewVersion(NoWCEGAL10num0[:,withinrange])
NoWCEGAL10num1average = SumPerSiteNewVersion(NoWCEGAL10num1[:,withinrange])
NoWCEGAL10num2average = SumPerSiteNewVersion(NoWCEGAL10num2[:,withinrange])

fig = plt.figure(figsize = (14,14) ,dpi=200)

ax1 = fig.add_subplot(3,1,1)
ax1.step(withinrange, NoWCEGAL10num0average, label = 'Cluster 1', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax2 = fig.add_subplot(3,1,2)
ax2.step(withinrange, NoWCEGAL10num1average, label = 'Cluster 2', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax3 = fig.add_subplot(3,1,3)
ax3.step(withinrange, NoWCEGAL10num2average, label = 'Cluster 3', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.xlabel('i')
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

fig.savefig(outputdir + 'AverageNoWCEGAL10.pdf', dpi=1200, transparent=True)
plt.show()

# figure GAL WCE average clusters for GAL10
GALWCEGAL10num0average = SumPerSiteNewVersion(GALWCEGAL10num0[:,withinrange])
GALWCEGAL10num1average = SumPerSiteNewVersion(GALWCEGAL10num1[:,withinrange])
GALWCEGAL10num2average = SumPerSiteNewVersion(GALWCEGAL10num2[:,withinrange])

fig = plt.figure(figsize = (14,14) ,dpi=200)

ax1 = fig.add_subplot(3,1,1)
ax1.step(withinrange, GALWCEGAL10num0average, label = 'Cluster 1', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax2 = fig.add_subplot(3,1,2)
ax2.step(withinrange, GALWCEGAL10num1average, label = 'Cluster 2', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax3 = fig.add_subplot(3,1,3)
ax3.step(withinrange, GALWCEGAL10num2average, label = 'Cluster 3', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.xlabel('i')
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

fig.savefig(outputdir + 'AverageGALWCEGAL10.pdf', dpi=1200, transparent=True)
plt.show()

# figure RAF WCE average clusters for GAL10
RAFWCEGAL10num0average = SumPerSiteNewVersion(RAFWCEGAL10num0[:,withinrange])
RAFWCEGAL10num1average = SumPerSiteNewVersion(RAFWCEGAL10num1[:,withinrange])
RAFWCEGAL10num2average = SumPerSiteNewVersion(RAFWCEGAL10num2[:,withinrange])

fig = plt.figure(figsize = (14,14) ,dpi=200)

ax1 = fig.add_subplot(3,1,1)
ax1.step(withinrange, RAFWCEGAL10num0average, label = 'Cluster 1', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax2 = fig.add_subplot(3,1,2)
ax2.step(withinrange, RAFWCEGAL10num1average, label = 'Cluster 2', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax3 = fig.add_subplot(3,1,3)
ax3.step(withinrange, RAFWCEGAL10num2average, label = 'Cluster 3', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.xlabel('i')
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

fig.savefig(outputdir + 'AverageRAFWCEGAL10.pdf', dpi=1200, transparent=True)
plt.show()


# In[206]:


motif = seqGal1
plotrange, withinrange = llr_in_range(new_indexplus, fasta, motif)

# figure No WCE average clusters for GAL10
NoWCEGAL1num0average = SumPerSiteNewVersion(NoWCEGAL1num0[:,withinrange])
NoWCEGAL1num1average = SumPerSiteNewVersion(NoWCEGAL1num1[:,withinrange])
NoWCEGAL1num2average = SumPerSiteNewVersion(NoWCEGAL1num2[:,withinrange])
NoWCEGAL1num3average = SumPerSiteNewVersion(NoWCEGAL1num3[:,withinrange])

fig = plt.figure(figsize = (14,14) ,dpi=200)

ax1 = fig.add_subplot(4,1,1)
ax1.step(withinrange, NoWCEGAL1num0average, label = 'Cluster 1', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.xlabel('i')
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax2 = fig.add_subplot(4,1,2)
ax2.step(withinrange, NoWCEGAL1num1average, label = 'Cluster 2', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax3 = fig.add_subplot(4,1,3)
ax3.step(withinrange, NoWCEGAL1num2average, label = 'Cluster 3', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax4 = fig.add_subplot(4,1,4)
ax4.step(withinrange, NoWCEGAL1num3average, label = 'Cluster 4', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.xlabel('i')
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

fig.savefig(outputdir + 'AverageNoWCEGAL1.pdf', dpi=1200, transparent=True)
plt.show()

# figure GAL WCE average clusters for GAL10
GALWCEGAL1num0average = SumPerSiteNewVersion(GALWCEGAL1num0[:,withinrange])
GALWCEGAL1num1average = SumPerSiteNewVersion(GALWCEGAL1num1[:,withinrange])
GALWCEGAL1num2average = SumPerSiteNewVersion(GALWCEGAL1num2[:,withinrange])

fig = plt.figure(figsize = (14,14) ,dpi=200)

ax1 = fig.add_subplot(3,1,1)
ax1.step(withinrange, GALWCEGAL1num0average, label = 'Cluster 1', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax2 = fig.add_subplot(3,1,2)
ax2.step(withinrange, GALWCEGAL1num1average, label = 'Cluster 2', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax3 = fig.add_subplot(3,1,3)
ax3.step(withinrange, GALWCEGAL1num2average, label = 'Cluster 3', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.xlabel('i')
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

fig.savefig(outputdir + 'AverageGALWCEGAL1.pdf', dpi=1200, transparent=True)
plt.show()

# figure RAF WCE average clusters for GAL10
RAFWCEGAL1num0average = SumPerSiteNewVersion(RAFWCEGAL1num0[:,withinrange])
RAFWCEGAL1num1average = SumPerSiteNewVersion(RAFWCEGAL1num1[:,withinrange])
RAFWCEGAL1num2average = SumPerSiteNewVersion(RAFWCEGAL1num2[:,withinrange])
RAFWCEGAL1num3average = SumPerSiteNewVersion(RAFWCEGAL1num3[:,withinrange])

fig = plt.figure(figsize = (14,14) ,dpi=200)

ax1 = fig.add_subplot(4,1,1)
ax1.step(withinrange, RAFWCEGAL1num0average, label = 'Cluster 1', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax2 = fig.add_subplot(4,1,2)
ax2.step(withinrange, RAFWCEGAL1num1average, label = 'Cluster 2', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax3 = fig.add_subplot(4,1,3)
ax3.step(withinrange, RAFWCEGAL1num2average, label = 'Cluster 3', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

ax4 = fig.add_subplot(4,1,4)
ax4.step(withinrange, RAFWCEGAL1num3average, label = 'Cluster 4', linewidth = 1, linestyle = '-') # scatter, s = 50, marker = 'o')
plt.ylim((-1,1))
plt.xlabel('i')
plt.legend(loc='upper right')
plt.grid(True)
plt.ylabel('average LLR')
plt.xlim(plotrange[0],plotrange[1])

fig.savefig(outputdir + 'AverageRAFWCEGAL1.pdf', dpi=1200, transparent=True)
plt.show()


# In[ ]:




