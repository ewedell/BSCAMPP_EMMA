#from dendropy import *
import numpy as np
import heapq
import treeswift
import itertools
from collections import deque
from os.path import expanduser,isfile
import random
import statistics
import copy

# store bracket open/close for convenience in label parsing
BRACKET = {
    '[': ']', # square bracket
    '{': '}', # curly bracket
    "'": "'", # single-quote
    '"': '"', # double-quote
}


#separete the query and ref sequence from the alignment file

def read_data(aln):
    """ Load the query and reference sequence from the alignment file
    
    Parameters
    ----------
    aln : multiple sequence alignment containing reference taxa and query sequence
    
    Returns
    -------
    dictionary containing sequences with taxon label keys
    
    """

    f = open(aln)
    result = dict()

    taxa = ""
    seq = ""
    for line in f:
        if line[0] == '>':
            if taxa != "":
                result[taxa] = seq
            taxa = line[1:-1]
            seq = ""

        elif line == "/n":
            continue
        else:
            seq += line[:-1]

    if taxa != "":
        result[taxa] = seq


    return result

def split_alignment(aln_dict, leaf_dict):
    """ Separate the query sequences from the reference sequences
    
    Parameters
    ----------
    aln_dict : Sequence dictionary with taxon label keys
    leaf_dict : Sequence dictionary with leaf label keys (queries are not in backbone tree)
    
    Returns
    -------
    separate dictionaries containing query sequences and referece sequences with taxon label keys
    
    """
    ref = dict()
    query = dict()

    for key, value in aln_dict.items():
        if key not in leaf_dict:
            query[key] = value
        else:
            ref[key] = value

    return ref, query
    

def write_fasta(aln_path, aln_dict, aligned=True):
        
        f = open(aln_path, "w")
        for label, seq in aln_dict.items():
            if label != '':
                f.write(">"+label+"\n")
                if aligned:
                    f.write(seq+"\n")
                else:
                    f.write(seq.replace('-','')+"\n")
        f.close()


def subtree_nodes_with_edge_length(tree, leaf_y, n):
    """ Returns list of length n of leaves closest to sister taxon (minimizing edge weights)
    
    Parameters
    ----------
    tree : treeswift tree object
    leaf_y : treeswift node for closest sister taxon
    n = number of taxa contained in subtree
    
    Returns
    -------
    list of taxon labels corresponding to leaves in the subtree
    """
    queue = [(leaf_y.get_edge_length(), leaf_y.get_parent())]

    leaves = [leaf_y]
    visited = {leaf_y}

    while len(leaves) < n:
        try:
            (length, node) = heapq.heappop(queue)
        except IndexError:
            break

        visited.add(node)
        if node.is_leaf() and node.get_label() != '':
            leaves.append(node)

        adjacent_nodes = node.child_nodes()
        if not node.is_root():
            adjacent_nodes.append(node.get_parent())

        for neighbor in adjacent_nodes:
            if neighbor not in visited:
                if neighbor == node.get_parent():
                    heapq.heappush(queue, (length+node.get_edge_length(), neighbor))
                else:
                    heapq.heappush(queue, (length+neighbor.get_edge_length(), neighbor))

    result = []
    for item in leaves:
        result.append(item.get_label())

    return result

