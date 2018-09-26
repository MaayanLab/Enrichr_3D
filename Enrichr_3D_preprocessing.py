import os
import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt
import plotly.plotly as py

import pickle

from sklearn.metrics.pairwise import pairwise_distances
from sklearn.manifold import TSNE



from gensim import corpora, models, similarities
from gensim.corpora import SvmLightCorpus



def readGMT(path):
    with open(path) as file:
        FILE = file.readlines()
        FILE  = [[x.split(",")[0] for x in line.split("\t") if x != '\n'] for line in FILE]
    return FILE

def GMT_dict():
    gmt_dict = {}
    try:
        os.remove('libs/.DS_Store')
    except:
        pass
    for file in os.listdir("libs"):
        lib_name = file.strip(".txt")
        gmt_dict[lib_name] = readGMT("libs/"+file)
    return gmt_dict
gmt_dict = GMT_dict()




def get_total_genes(gmt_dict):
    uniqueGenes = set()
    for lib in gmt_dict:
        for sig in gmt_dict[lib]:
            uniqueGenes = set(sig[2:]).union(uniqueGenes)
    # genes = [[set(sig[2:]) for sig in gmt_dict[lib]] for lib in gmt_dict]
    print("There are %d unique genes across all libraries." % len(uniqueGenes))
    return sorted(list(uniqueGenes))
geneUniverse = get_total_genes(gmt_dict)



def GMTdict_to_binaryDF(gmt_dict, geneUniverse):
    def assign_binary(gene_fromUniverse):
        bit = 1 if gene_fromUniverse in genes else 0
        return bit
    binary_dict = {}
    for lib in gmt_dict:
        print("Processing library : %s" % lib)
        for sig in gmt_dict[lib]:
            compound_name = lib+" @ "+sig[0]
            genes = list(set([gene.upper() for gene in sig[2:]]))
            binary_dict[compound_name] = list(map(assign_binary, geneUniverse))
    return pd.DataFrame(binary_dict)
binary_df = GMTdict_to_binaryDF(gmt_dict, geneUniverse)
binary_df.head()
binary_df.sum() # Get number of genes per signature
print("Average number of genes per signature = %d" % (binary_df.sum().sum() / binary_df.sum().__len__()))
binary_df.to_csv('Enrichr_binary_matrix.csv', index=False)




##### Jaccard Similarity Index #####
jac_sim = 1- pairwise_distances(binary_df.T, metric = "hamming")
print(jac_sim.shape)

jac_df = pd.DataFrame(jac_sim, index=binary_df.columns, columns=binary_df.columns)
jac_df.head()
print("Min distance = %.2f" % jac_df.min().min())
print("Max distance = %.2f" % jac_df.max().max())
jac_df.to_csv("Enrichr_JaccardIndex_normalized.csv")

# Normalize to make 0-1
jac_norm_df = (jac_df - jac_df.mean()) / (jac_df.max() - jac_df.min())
print("Min distance = %.2f" % jac_norm_df.min().min())
print("Max distance = %.2f" % jac_norm_df.max().max())
jac_norm_df.to_csv("Enrichr_JaccardIndex_normalized.csv")


# Plot Jaccard Index (2D)
sn.clustermap(jac_df)
plt.show()



#### Dimensionality Reduction ####

#### ---------- Latent Semantic Indexing (LSI) ---------- ####
# LSI as PCA for binary or multinomial data
## Convert Data Using : Term Frequency - Inverse Document Frequency
# with open("Enrichr_binary_matrix.csv") as file:
#     binary_txt = file.readlines()
#
# def GMTdict_to_corpus(gmt_dict):
#     corpus = []; corpus_names =[];
#     for lib in gmt_dict:
#         for sig in gmt_dict[lib]:
#             sig_entry = tuple(sorted(sig[2:]))
#             corpus.append(sig_entry)
#             corpus_names.append(lib+" @ "+sig[0])
#
# import gensim.downloader as api
# dataset = api.load("text8")
#
# dct = Dictionary(binary_txt)
# corpus = SvmLightCorpus(binary_txt)
# tfidf = models.TfidfModel( corpus )
# c_tfidf = tfidf[corpus]
#
# lsi = models.LsiModel( c_tfidf, id2word = None, num_topics = num_topics )
#


#### ---------- t-SNE ---------- ####

def run_tSNE (df, n_components, save=True):
    tsne_results = TSNE(n_components=n_components, verbose=True).fit(df)
    tsne_data = tsne_results.fit_transform(df)
    if save==True:
        pickle.dump( [tsne_results, tsne_data], open( "tSNE_results.p", "wb" ) )
    return tsne_results, tsne_data
binary_df_T = binary_df.transpose()
tsne_results, tsne_data = run_tSNE(binary_df_T, 3, True)

tsne_df = pd.DataFrame({'tsne1':tsne_data[:,0], 'tsne2':tsne_data[:,1], 'tsne3':tsne_data[:,2]}, index= binary_df_T.index)
tsne_df['lib_names'] = tsne_df.index.str.split(" @ ").str[0].tolist()

tsne_pickle = pickle.load( open( "tSNE_results.p", "rb" ) )

# 2D tSNE plot
plt.subplots()
plt.title('t-SNE\nBinary Gene Signatures from Enrichr')
plt.xlabel('Dimension 1')
plt.ylabel('Dimension 2')
plt.scatter(x=tsne_df['tsne1'], y=tsne_df['tsne2'], alpha=0.5)
plt.show()


# 3D tSNE plot
def tSNE_3D(df):
    data = []
    # clusters = []
    nGroups = len(df['lib_names'].unique())
    import colorlover as cl
    colors = cl.scales[str(nGroups)]['qual']['Paired']
    # colors = ['rgb(228,26,28)', 'rgb(55,126,184)', 'rgb(77,175,74)']


    for i in range(nGroups):
        name = df['lib_names'].unique()[i]
        color = colors[i]
        x = df[df['lib_names'] == name]['tsne1']
        y = df[df['lib_names'] == name]['tsne2']
        z = df[df['lib_names'] == name]['tsne3']

        trace = dict(
            name=name,
            x=x, y=y, z=z,
            type="scatter3d",
            mode='markers',
            marker=dict(size=3, color=color, line=dict(width=0)))
        data.append(trace)

    layout = dict(
        width=800,
        height=550,
        autosize=False,
        title='t-SNE\nBinary Gene Signatures from Enrichr',
        scene=dict(
            xaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            yaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            zaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            aspectratio=dict(x=1, y=1, z=0.7),
            aspectmode='manual'
        ),
    )

    fig = dict(data=data, layout=layout)
    url = py.plot(fig, filename='Enrichr_3D_test', validate=False)
    return url

url = tSNE_3D(tsne_df)