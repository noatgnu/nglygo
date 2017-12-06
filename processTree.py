import argparse
import tempfile
from Bio import Phylo

parser = argparse.ArgumentParser()
parser.add_argument("-ct", type=str, dest="ct")
parser.add_argument("-asrt", type=str, dest="asrt")
parser.add_argument("-o", type=str, dest="out")
a = parser.parse_args()
currentTree = Phylo.read(a.ct, format="newick")
reconstructed = Phylo.read(a.asrt, format="newick")
for c1, c2 in zip(currentTree.find_clades(), reconstructed.find_clades()):
    if not c1.name:
        c1.name = str(c2.confidence)
Phylo.write(currentTree, a.out, 'phyloxml')
with tempfile.NamedTemporaryFile(mode='w+t') as temp:
    Phylo.write(currentTree, temp, 'newick')
    temp.flush()
