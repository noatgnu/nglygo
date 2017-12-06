import argparse

from Bio import Phylo

parser = argparse.ArgumentParser()
parser.add_argument("-ct", type=str, dest="ct")
parser.add_argument("-asrt", type=str, dest="asrt")
parser.add_argument("-o", type=str, dest="out")

print("test")
a = parser.parse_args()
currentTree = Phylo.read(a.ct, format="newick")
reconstructed = Phylo.read(a.asrt, format="newick")
for c1, c2 in zip(currentTree.find_clades(), reconstructed.find_clades()):
    print(c1.name,c2.name)

with open(a.out, "wt", newline="") as outFile:
    Phylo.write(currentTree, outFile, "newick")
    outFile.flush()
