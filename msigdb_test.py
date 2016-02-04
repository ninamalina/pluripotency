from collections import defaultdict
from os.path import exists
import pickle
import networkx as nx
import numpy as np
from Data.get_names import get_dict
import os
import random
import time
from sklearn import metrics

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


# tested = ["ECAT1", "DPPA5", "FBXO15", "NANOG", "ERAS", "DNMT3l", "ECAT8", "GDF3", "SOX15", "DPPA4", "DPPA2", "FTHL17",
#               "SALL4", "REX1", "UTF1", "TCL1", "DPPA3",  "B-CATENIN",  "STAT3", "GRB2"]
# simon = ["CD4", "IL2RA", "PTPRC", "FOXP3", "CTLA4", "TIGIT", "CD40LG", "THEMIS", "SATB1"]
#
# four = ["OCT3", "OCT4", "OCT3/4", "SOX2","KLF4", "C-MYC"]
#
# genes_with_gos = pickle.load(open('names_to_test.p', 'rb'))



time0 = time.clock()
interactor_name_dict = {}
interactor_index_dict = {}

# check if Pajek NET file already exists
if exists(os.path.expanduser("~/PycharmProjects/pluripotency/Data/Homo_sapiens_interactions.net")):
#if exists( r'.\Data\Homo_sapiens_interactions.net'):
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
if exists("gene_neighbours.p"):
    gene_neighbours = pickle.load(open("gene_neighbours.p",'rb'))
else:
    for gene in name_interactor_dict:
        # print(gene, n_neighbor(G, name_interactor_dict[temp[gene]],2))
        gene_neighbours[name_interactor_dict[gene]] = nx.neighbors(G,name_interactor_dict[gene])

    pickle.dump(gene_neighbours, open('gene_neighbours.p', "wb"))

print("gene_neighbours computed")

gene_neighbours_2 = defaultdict(list)
if exists("gene_neighbours_2.p"):
    gene_neighbours_2 = pickle.load(open("gene_neighbours_2.p",'rb'))
else:
    for gene in name_interactor_dict:
        # print(gene, n_neighbor(G, name_interactor_dict[temp[gene]],2))
        for gen1 in gene_neighbours[name_interactor_dict[gene]]:
            gene_neighbours_2[name_interactor_dict[gene]] += nx.neighbors(G,gen1)
        # gene_neighbours[name_interactor_dict[temp[gene]]] = nx.neighbors(G,name_interactor_dict[temp[gene]])
    for gene in gene_neighbours_2:
        gene_neighbours_2[gene] = set(gene_neighbours_2[gene])

    pickle.dump(gene_neighbours_2, open('gene_neighbours_2.p', "wb"))

print("gene_neighbours_2 computed")


l = 1000
n = 4

four_random_scores= []

for i in range(l):
    four_random = random.sample(list(name_interactor_dict), n)
    s = set(gene_neighbours_2[name_interactor_dict[four_random[0]]])
    u = set(gene_neighbours_2[name_interactor_dict[four_random[0]]])
    for a in four_random:
        s = s.intersection(set(gene_neighbours_2[name_interactor_dict[a]]))
        u = u.union(set(gene_neighbours_2[name_interactor_dict[a]]))
    if len(u) > 0 :
        m = len(s) / float(len(u))
        four_random_scores.append(m)


V = G.nodes()
degrees  = sorted([G.degree(node) for node in V])


if len(degrees) % 2 == 0:
    d_m = (degrees[int(len(degrees) / 2)] + degrees[int(len(degrees) / 2 - 1)]) / 2.
else:
    d_m = degrees[int(len(degrees) / 2)]

print("d_m:", d_m)

path = 'Data/msigdb_genesets/'
listing = os.listdir(path)

profiles_file = open("profiles.txt", "w+")
profiles_file.write("n_S, m_S, n, m, c_S, internal_density, edges_inside, average_degree, fomd, tpr, expansion, cut_ratio, conductance, normalized_cut, maximum_odf, average_odf, flake_odf, modularity, auc_score\n")

k = 0

for infile in listing:
    f = open(os.path.expanduser("~/PycharmProjects/pluripotency/Data/msigdb_genesets" + "//" + infile))

    print("current file is: " + infile)

    msigdb_genes = []
    for line in f:
        splited = line.split("\t")[2:]
        if len(splited) > 5:
            msigdb_genes.append([a.strip() for a in splited])

    print(msigdb_genes)

    l = 1000
    n1 = 4
#moj graf plus njihovi sosedi


    geneset_scores = []

    for msigdb_gene_set in msigdb_genes:
        print(k)
        k += 1
        temp = get_dict(msigdb_gene_set) # name_officialName_dict
        if (len(temp) > 5 and len(temp) < 20) :

            s = list(set([name_interactor_dict[name] for name in list(temp) if name in name_interactor_dict.keys()]))
            s_with_neighbours = list(set([item for sublist in [gene_neighbours[interactor] for interactor in s if interactor in gene_neighbours.keys()] for item in sublist])) + s

            S_graph = G.subgraph(s)
            S_plus_graph = G.subgraph(s_with_neighbours)

            n_S = float(S_graph.number_of_nodes())
            m_S = float(S_graph.number_of_edges())
            n = float(G.number_of_nodes())
            m = float(G.number_of_edges())

            c_S = 0
            for (u,v) in S_plus_graph.edges():
                if (u in S_graph.nodes() and v not in S_graph.nodes()) or (u not in S_graph.nodes() and v in S_graph.nodes()):
                   c_S += 1

            # if (c_S == 0):
            #     print(S_graph.number_of_edges(), S_plus_graph.number_of_edges())
            #     print(S_graph.nodes())
            #     print("S plus nodes", S_plus_graph.nodes())
            #     print("S plus edges", S_plus_graph.edges())
            #     print("*", [(interactor, gene_neighbours[interactor]) for interactor in s if interactor in gene_neighbours.keys()])
            #     # c_S = 0
            #     # for (u,v) in S_plus_graph.edges():
            #     #     print(u, v)
            #     #     if (u in S_graph.nodes() and v not in S_graph.nodes()) or (u not in S_graph.nodes() and v in S_graph.nodes()):
            #     #         c_S += 1


            internal_density = 2 * m_S / (n_S * (n_S - 1))

            edges_inside = m_S

            average_degree = 2 * m_S / n_S


            fomd = sum(np.array([1 for node in S_graph.nodes() if S_graph.degree(node) > d_m])) / n_S

            tpr = 0

            for u in S_graph.nodes():
                u_count = 0
                for (v, w) in S_graph.edges():
                    if (u,v) in S_graph.edges() and (u,w) in S_graph.edges():
                        u_count += 1
                if u_count > 0:
                    tpr += 1
            tpr = tpr / n_S

            expansion = c_S / n_S

            cut_ratio = c_S / (n_S * (n - n_S))

            if (2 * m_S + c_S) > 0 :
                conductance = c_S / (2 * m_S + c_S)
            else:
                conductance = 0.0

            normalized_cut = conductance + c_S / (2 * (m - m_S) + c_S)


            odfs = [sum(np.array([1 for edge in S_plus_graph.edges(node) if edge not in S_graph.edges(node)])) / S_plus_graph.degree(node) for node in S_graph.nodes() if S_plus_graph.degree(node) !={} and S_plus_graph.degree(node) != 0]


            maximum_odf = 0.0
            average_odf = 0.0

            if len(odfs) > 0 :
                maximum_odf = max(odfs)
                average_odf = np.average(np.array(odfs))

            flake_odf = sum([1 for node in S_graph.nodes() if S_plus_graph.degree(node) != {} and S_graph.number_of_edges(node) < S_plus_graph.degree(node) / 2]) / n_S


            # modularity =

            degrees = []
            for node in S_graph.nodes():
                degree = 0
                if S_plus_graph.degree(node) != {}:
                    degree = S_plus_graph.degree(node)
                degrees.append(degree)



            # print("n_S:", n_S)
            # print("m_S:", m_S)
            # print("n:", n)
            # print("m:", m)
            # print("c_S:", c_S)
            # print("internal_density:", internal_density)
            # print("edges_inside:", edges_inside)
            # print("average_degree:", average_degree)
            # print("fomd:", fomd)
            # print("tpr:", tpr)
            # print("expansion:", expansion)
            # print("cut_ratio:", cut_ratio)
            # print("conductance:", conductance)
            # print("normalized_cut:", normalized_cut)
            # print("maximum_odf:", maximum_odf)
            # print("average_odf:", average_odf)
            # print("flake_odf:", flake_odf)



            #upostevas tud sosede ampak une povezave k majo sosedi vn, ne upostevas
            degrees_neighbours = []

            for u in S_plus_graph.nodes():
                if u not in S_graph.nodes():
                    d = 0
                    for (u,v) in S_plus_graph.edges(u):
                        if v in S_graph.nodes() or u in S_graph.nodes():
                            d += 1
                    if d > 0:
                        degrees_neighbours.append(d)
            # print(sum(degrees_neighbours))
            # print(degrees)
            # print(degrees_neighbours)
            a = []
            for i in range(100):
                rg = nx.configuration_model(degrees + degrees_neighbours, create_using=None, seed=None)
                a.append(rg.subgraph(range(len(degrees))).number_of_edges())
            e_m = sum(a) / 100.

            modularity = (m_S - e_m) / 4.


            four_real_scores = []

            for i in range(l):
                four_real = random.sample(list(temp), n1)

                comb_num_dict = {}
                s = set(gene_neighbours_2[name_interactor_dict[temp[four_real[0]]]])
                u = set(gene_neighbours_2[name_interactor_dict[temp[four_real[0]]]])

                for a in four_real:
                    s = s.intersection(set(gene_neighbours_2[name_interactor_dict[temp[a]]]))
                    u = u.union(set(gene_neighbours_2[name_interactor_dict[temp[a]]]))

                # if(c_S == 0):
                #     print(s, u)

                if len(u) > 0:
                    jc_score = len(s) / float(len(u))
                    four_real_scores.append(jc_score)


            if(c_S == 0):
                print("c_s = 0, four_real_scores = ", four_real_scores)


            y = np.array([0]*len(four_random_scores) + [1]*len(four_real_scores))
            pred = np.array(four_random_scores + four_real_scores)
            fpr, tspr, thresholds = metrics.roc_curve(y, pred, pos_label=1)
            auc_score = metrics.auc(fpr, tspr)
            geneset_scores.append(auc_score)

            profile = [str(x) for x in [n_S, m_S, n, m, c_S, internal_density, edges_inside, average_degree, fomd, tpr, expansion, cut_ratio, conductance, normalized_cut, maximum_odf, average_odf, flake_odf, modularity, auc_score]]
            s = ",".join(profile) + "\n"
            profiles_file.write(s)
profiles_file.close()
print("time in minutes: ", (time.clock() - time0) / 60)
    #
    # print(geneset_scores)
    # print(len(geneset_scores))
    # print(np.array(geneset_scores).mean())









# plt.figure()
# plt.plot(range(len(scores)), sorted(scores), 'ro')
# plt.savefig("bf_scores.pdf")
# plt.close()


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

