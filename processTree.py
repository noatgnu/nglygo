import argparse

from Bio import Phylo

parser = argparse.ArgumentParser()
parser.add_argument("-ct", type=str, dest="ct")
parser.add_argument("-asrt", type=str, dest="asrt")
parser.add_argument("-o", type=str, dest="out")
if "__name__" == "__main__":
    a = parser.parse_args()
    currentTree = Phylo.read(a.ct)
    reconstructed = Phylo.read(a.asrt)
    for c1, c2 in zip(currentTree.find_clades(order="level"), reconstructed.find_clades(order="level")):
        if "_" not in c2.name:
            c1.name = c2.name
    with open(a.out, "wt", newline="") as outFile:
        Phylo.write(currentTree, outFile, "newick")
        outFile.flush()
