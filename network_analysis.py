from collections import defaultdict
from os.path import exists
import pickle
import networkx as nx
import matplotlib.pylab as plt
import numpy as np
from sklearn import manifold
from Data.get_names import get_dict
import os
from itertools import combinations
import operator

__author__ = 'NinaHP'

def distance_matrix(genes): # pairwise jaccard similarity

    print(genes)
    print(len(genes))

    dm = np.zeros((len(genes), len(genes)))
    for i in range(len(genes)):
        for j in range(len(genes)):
            presek = set(gene_neighbours_2[name_interactor_dict[temp[genes[i]]]]).intersection(
                set(gene_neighbours_2[name_interactor_dict[temp[genes[j]]]]))
            unija = set(gene_neighbours_2[name_interactor_dict[temp[genes[i]]]]).union(
                set(gene_neighbours_2[name_interactor_dict[temp[genes[j]]]]))
            dm[i][j] = 1 - len(presek) / float(len(unija))

    return dm


tested = ["ECAT1", "DPPA5", "FBXO15", "NANOG", "ERAS", "DNMT3l", "ECAT8", "GDF3", "SOX15", "DPPA4", "DPPA2", "FTHL17",
              "SALL4", "REX1", "UTF1", "TCL1", "DPPA3",  "B-CATENIN",  "STAT3", "GRB2"]
simon = ["CD4", "IL2RA", "PTPRC", "FOXP3", "CTLA4", "TIGIT", "CD40LG", "THEMIS", "SATB1"]

four = ["OCT3", "OCT4", "OCT3/4", "SOX2","KLF4", "C-MYC"]

genes_with_gos = pickle.load(open('names_to_test.p', 'rb'))

temp = get_dict(tested + simon + four + genes_with_gos) # name_officialName_dict

# temp = {'OCT3':'POU5F1', 'SOX2':'SOX2', 'KLF4':'KLF4', 'C-MYC':'MYC', 'DPPA4':'DPPA4', 'FBXO15':'FBXO15', 'NANOG':'NANOG', 'ERAS':'ERAS', 'DNMT3l': 'DNMT3L', 'GDF3':'GDF3', 'SOX15':'SOX15', 'DPPA2':'DPPA2', 'SALL4':'SALL4', 'REX1':'ZFP42', 'UTF1':'UTF1', 'TCL1':'TCL1A', 'STAT3':'STAT3', 'GRB2':'GRB2', 'CD4':'CD4', 'IL2RA':'IL2RA', 'PTPRC':'PTPRC', 'FOXP3':'FOXP3', 'CTLA4':'CTLA4', 'TIGIT':'TIGIT', 'CD40LG':'CD40LG', 'THEMIS':'THEMIS', 'SATB1':'SATB1'}
# temp = {'OCT3':'POU5F1', 'SOX2':'SOX2', 'KLF4':'KLF4', 'C-MYC':'MYC', 'DPPA4':'DPPA4', 'DNMT3l': 'DNMT3L', 'UTF1':'UTF1', 'TCL1':'TCL1A', 'STAT3':'STAT3', 'GRB2':'GRB2', 'CD4':'CD4'}

plt.savefig('labels.png')

interactor_name_dict = {}
interactor_index_dict = {}

# check if Pajek NET file already exists
if exists(os.path.expanduser("~/PycharmProjects/pluripotency/Data/Homo_sapiens_interactions.net")):
#if exists( r'.\Data\Homo_sapiens_interactions.net'):
    print("tuki")
    interactor_name_dict = pickle.load( open(os.path.expanduser("~/PycharmProjects/pluripotency/Data/interactor_name_dict.p"), "rb" ) )
    interactor_index_dict = pickle.load( open(os.path.expanduser("~/PycharmProjects/pluripotency/Data/interactor_index_dict.p"), "rb" ) )

    # interactor_name_dict = pickle.load( open( r'.\Data\interactor_name_dict.p', "rb" ) )
    # interactor_index_dict = pickle.load( open( r'.\Data\interactor_index_dict.p', "rb" ) )

else:
    # read interactions
    #f = open(r'.\Data\BIOGRID-ORGANISM-Homo_sapiens-3.4.129.tab.txt', 'r')
    f = open(os.path.expanduser("~/PycharmProjects/pluripotency/Data/BIOGRID-ORGANISM-Homo_sapiens-3.4.129.tab.txt"))
    interactions = []
    # za nove podatke je 36 !!
    # za stare podatke je 32 !!
    for _ in range(36):
        next(f)
    for line in f:
        line = line.split('\t')
        first_interactor = line[0].strip()
        second_interactor = line[1].strip()
        if first_interactor not in interactor_name_dict:
            interactor_name_dict[first_interactor] = line[2].strip()
        if second_interactor not in interactor_name_dict:
            interactor_name_dict[second_interactor] = line[3].strip()
        interactions.append((first_interactor, second_interactor))

    f.close()
    # pickle.dump(interactor_name_dict, open( r'.\Data\interactor_name_dict.p', "wb" ) )
    pickle.dump(interactor_name_dict, open(os.path.expanduser("~/PycharmProjects/pluripotency/Data/interactor_name_dict.p"), "wb" ) )

    # create pajek NET file
    f = open(os.path.expanduser("~/PycharmProjects/pluripotency/Data/Homo_sapiens_interactions.net"), 'w+')
    # f = open(r'.\Data\Homo_sapiens_interactions.net', 'w+')

    f.write('*Vertices '+str(len(interactor_name_dict))+'\n')
    i = 1
    for interactor in interactor_name_dict.keys():
        f.write(str(i)+' "'+ interactor+ '"\n')
        interactor_index_dict[interactor] = i
        i += 1

    # pickle.dump(interactor_index_dict, open( r'.\Data\interactor_index_dict.p', "wb" ) )
    pickle.dump(interactor_index_dict, open(os.path.expanduser("~/PycharmProjects/pluripotency/Data/interactor_name_dict.p"), "wb" ) )

    f.write('*Edges\n')
    for interaction in interactions:
        f.write(str(interactor_index_dict[interaction[0]])+' '+str(interactor_index_dict[interaction[1]])+'\n')

    f.close()



# G = nx.read_pajek(r'.\Data\Homo_sapiens_interactions.net')
G = nx.read_pajek(os.path.expanduser("~/PycharmProjects/pluripotency/Data/Homo_sapiens_interactions.net"))

# nx.write_edgelist(G, r'.\Data\edge_list.txt')

# betweenes_centrality = nx.betweenness_centrality(G)
# pickle.dump(betweenes_centrality, open(r'.\Data\betweenes_centrality_dict_2012.p', 'wb'))
# degree_centrality = nx.degree_centrality(G)
# pickle.dump(degree_centrality, open(r'.\Data\degree_centrality_dict_2012.p', 'wb'))
# closenes_centrality = nx.closeness_centrality(G)
# pickle.dump(closenes_centrality, open(r'.\Data\closenes_centrality_dict_2012.p', 'wb'))
# page_rank_centrality = nx.pagerank_numpy(G)


# betweenes_centrality = pickle.load(open(r'.\Data\betweenes_centrality_dict.p', 'rb'))
# degree_centrality = pickle.load(open(r'.\Data\degree_centrality_dict.p','rb'))
# closenes_centrality = pickle.load(open(r'.\Data\closenes_centrality_dict.p', 'rb'))

betweenes_centrality = pickle.load(open(os.path.expanduser("~/PycharmProjects/pluripotency/Data/betweenes_centrality_dict.p"), 'rb'))
degree_centrality = pickle.load(open(os.path.expanduser("~/PycharmProjects/pluripotency/Data/degree_centrality_dict.p"),'rb'))
closenes_centrality = pickle.load(open(os.path.expanduser("~/PycharmProjects/pluripotency/Data/closenes_centrality_dict.p"), 'rb'))


name_interactor_dict = dict(zip(interactor_name_dict.values(),interactor_name_dict.keys()))
index_interactor_dict = dict(zip(interactor_index_dict.values(),interactor_index_dict.keys()))
tabela_genov = sorted(list(temp.keys()))


# shortest_paths = [list([-1]*len(temp)) for i in range(len(temp))]
#
# for i in range(len(temp)):
#     for j in range(len(temp)):
#         if nx.has_path(G,name_interactor_dict[temp[tabela_genov[i]]], name_interactor_dict[temp[tabela_genov[j]]]):
#             shortest_paths[i][j] = len(nx.shortest_path(G,name_interactor_dict[temp[tabela_genov[i]]], name_interactor_dict[temp[tabela_genov[j]]]))
#         else:
#             shortest_paths[i][j] = 9999
#         # print(len(nx.shortest_path(G,name_interactor_dict[temp[tabela_genov[i]]], name_interactor_dict[temp[tabela_genov[j]]])))


gene_neighbours = {}
for gene in tabela_genov:
    # print(gene, n_neighbor(G, name_interactor_dict[temp[gene]],2))
    gene_neighbours[name_interactor_dict[temp[gene]]] = nx.neighbors(G,name_interactor_dict[temp[gene]])

gene_neighbours_2 = defaultdict(list)
for gene in tabela_genov:
    # print(gene, n_neighbor(G, name_interactor_dict[temp[gene]],2))
    for gen1 in gene_neighbours[name_interactor_dict[temp[gene]]]:
        gene_neighbours_2[name_interactor_dict[temp[gene]]] += nx.neighbors(G,gen1)
    # gene_neighbours[name_interactor_dict[temp[gene]]] = nx.neighbors(G,name_interactor_dict[temp[gene]])
for gene in gene_neighbours_2:
    gene_neighbours_2[gene] = set(gene_neighbours_2[gene])


# print combinations of n genes, orderd by jaccard similarity
# for n in [3,4,5]:
#     comb_num_dict = {}
#     for c in combinations([key for key in temp if key  in simon],n):
#         s = set(gene_neighbours_2[name_interactor_dict[temp[c[0]]]])
#         u = set(gene_neighbours_2[name_interactor_dict[temp[c[0]]]])
#         for i in range(n):
#             s = s.intersection(set(gene_neighbours_2[name_interactor_dict[temp[c[i]]]]))
#             u = u.union(set(gene_neighbours_2[name_interactor_dict[temp[c[i]]]]))
#         comb_num_dict[c] = len(s) / float(len(u))
#
#     total_sort = sorted(comb_num_dict.items(), key=operator.itemgetter(1), reverse=True)
#     print('num of combination:', len(total_sort))
#     d = defaultdict(int)
#     for c in total_sort:
#         for g in c[0]:
#             d[g] += c[1]
#         print("*".join(c[0]), round(c[1],4))
#     print(d)


# dm = distance_matrix([key for key in temp if key in simon])

# plots genes,distance is jaccard similarity
# mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
# results = mds.fit(dm)
# coords = results.embedding_
# plt.subplots_adjust(bottom = 0.1)
# plt.scatter(
#     coords[:, 0], coords[:, 1], marker = 'o', label = [key for key in temp if key in simon]
#     )
#
# for label, x, y in zip([key for key in temp if key in simon], coords[:, 0], coords[:, 1]):
#     plt.annotate(
#         label,
#         xy = (x, y), xytext = (-1, 1),
#         textcoords = 'offset points', ha = 'right', va = 'bottom'
#         )
# plt.show()

# prints genes, orderd bz centrality measure
# a = {}
# for g in temp:
#     print(hits[1][name_interactor_dict[temp[g]]])
#     a[g] = hits[1][name_interactor_dict[temp[g]]]
# ranks = {}
# sort = sorted(a.items(), key=operator.itemgetter(1), reverse=True)
# for i in range(len(sort)):
#     ranks[sort[i][0]] = i+1
# for g in temp:
#     print(ranks[g])

print("Done!;)")