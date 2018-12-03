import argparse
from copy import deepcopy
import dendropy

p = argparse.ArgumentParser()
p.add_argument("-in", type=str, dest="i")
p.add_argument("-org", type=str, dest="org")
p.add_argument("-r", type=str, dest="r")
p.add_argument("-o", type=str, dest='o')
a = p.parse_args()


def tree_traverse(node, parent=None):
    for child in node.child_nodes():
        yield child
        yield from tree_traverse(child)


def decorate_tree(tree, replace_dict):
    t = deepcopy(tree)
    for n in tree_traverse(t.seed_node):
        if n.parent_node is not None:
            if n.label is None:
                if n.taxon != None:
                    if n.taxon._label in replace_dict:
                        tax = dendropy.Taxon(label=replace_dict[n.taxon._label])
                        n.taxon = tax
                        t.taxon_namespace.add_taxon(tax)
    return t


org = a.org.split(';')
replace_dict = dict()
if a.r is not '':
    its = a.r.split(';')
    for i in its:
        r = i.split(':')
        replace_dict[r[0]] = r[1]
dentree = dendropy.Tree.get_from_string(a.i, schema='newick')
retains = []
for n in tree_traverse(dentree.seed_node):
    if n.parent_node is not None:
        if n.label is None:
            if n.taxon != None:
                if n.taxon._label in org:
                    retains.append(n.taxon._label)
dentree.retain_taxa_with_labels(org)
new_tree = decorate_tree(dentree, replace_dict)
with open(a.o, 'w') as out_tree:
    new_tree.write(file=out_tree, schema='newick')
