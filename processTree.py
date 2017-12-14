import argparse

from Bio import Phylo

p = argparse.ArgumentParser()
p.add_argument("-ct", type=str, dest="ct")
p.add_argument("-asrt", type=str, dest="asrt")
p.add_argument("-b", type=str, dest="branches")
p.add_argument("-motifs", type=str, dest="motifs")
p.add_argument("-o", type=str, dest="out")
a = p.parse_args()
currentTree = Phylo.read(a.ct, format="newick")
reconstructed = Phylo.read(a.asrt, format="newick")

for c1, c2 in zip(currentTree.find_clades(), reconstructed.find_clades()):
    c1.test = "test"
    if not c1.name:
        c1.name = str(c2.confidence)
    else:
        c1.name = c2.name.split("_")[1]
Phylo.write(currentTree, a.out, 'newick')