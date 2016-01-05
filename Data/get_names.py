import os

__author__ = 'NinaHP'

def get_dict(list_of_names):
    d = {}
#    f = open(os.path.expanduser("~/Desktop/somefile.txt"))~/Desktop/somefile.t
    f = open(os.path.expanduser("~/PycharmProjects/pluripotency/Data/BIOGRID-ORGANISM-Homo_sapiens-3.4.129.tab.txt"))
    for _ in range(36):
        next(f)

    for line in f:
        line = line.split('\t')
        first_name = line[2].strip()
        second_name = line[3].strip()
        first_synonims = line[4].strip()
        second_synonims = line[5].strip()
        for a in first_synonims.split('|'):
            d[a.strip().upper()] = first_name
        for a in second_synonims.split('|'):
            d[a.strip().upper()] = second_name
        d[first_name.upper()] = first_name
        d[second_name.upper()] = second_name

    d_end = {}

    for g in list_of_names:
        if g in d:
            d_end[g] = d[g]
    return d_end