'''
New 
'''

import sys
import os
import utils
import shutil
import json
import time
import argparse
import treeswift
import copy
import merger
from collections import Counter

def main(args):
    tree_path = args.tree
    output_path = args.output_path
    temp_aln_path = args.temp_aln_path
    backbone_aln_path = args.backbone_aln_path
    n = args.subtreesize
    run = args.tmpfilenbr
    fragment_flag = args.fragmentflag
    queries_path = args.queries_path
    model = args.model
    max_query_size = args.max_query_size
    nbr_closest = args.votes

    # output path, ref, query, backbone tree, info
    t0 = time.perf_counter()
    tree = treeswift.read_tree_newick(tree_path)
    tree.resolve_polytomies()

    leaf_dict = tree.label_to_node(selection='leaves')
    print ('{} seconds loading tree'.format(time.perf_counter() - t0))

    # load up the alignments
    backbone_dict = utils.read_data(backbone_aln_path)
    
    if temp_aln_path == "":
        q_dict = utils.read_data(queries_path)
        # TO DO: add code to generate a temp_aln if not provided
    else:
        aln_dict = utils.read_data(temp_aln_path)
        ref_dict, q_dict = utils.split_alignment(aln_dict, leaf_dict)
    
    q_aln_path = "tmp{}/".format(run) + "qaln.fa"
    utils.write_fasta(q_aln_path, q_dict)       
        
    aln_path = "tmp{}/".format(run) + "aln.fa"
    utils.write_fasta(aln_path, ref_dict)    

    print ('{} seconds loading alignment'.format(time.perf_counter() - t0))

    # create 
    try:
        os.mkdir("tmp{}".format(run))
    except OSError as error:
    	pass
    #try:
    #    os.mkdir(output)
    #except OSError as error:
    #    pass


    # compute the hamming distance for each query sequence 
    # to each reference sequence to generate votes
    query_votes_dict = dict()
    query_top_vote_dict = dict()
    
    tmp_output_path = "tmp{}/".format(run) + "/closest.txt"

    if fragment_flag == True:
        os.system("./fragment_hamming {} {} {} {} {} {}".format(aln_path, len(ref_dict), q_aln_path, len(q_dict), tmp_output_path, nbr_closest))
    else:    
        os.system("./hamming {} {} {} {} {} {}".format(aln_path, len(ref_dict), q_aln_path, len(q_dict), tmp_output_path, nbr_closest))
    print ('{} seconds finding closest leaves'.format(time.perf_counter() - t0))
   
    f = open(tmp_output_path)
    for line in f:
        line = line.strip()
        y = line.split(',')
        name = y.pop(0)

        #print(name, y)
        for idx, taxon in enumerate(y):

            leaf, hamming = taxon.split(':')
            y[idx] = (leaf, int(hamming))

        y = sorted(y, key=lambda x: x[1])
        #print(y)
        for idx, taxon in enumerate(y):
            y[idx] = taxon[0]

        query_votes_dict[name] = y
        query_top_vote_dict[name] = y[0]
    f.close()

    print ('{} seconds processing closest leaves'.format(time.perf_counter() - t0))

    # tally votes
    lf_votes = Counter()
    leaf_queries = dict()
    for name, y in query_votes_dict.items():

        lf_votes.update(y)
        for ind, leaf in enumerate(y):
            top_vote = False
            if ind == 0:
                top_vote = True
            if leaf not in leaf_queries:           
                leaf_queries[leaf] = {(name,top_vote)}
            else:
                leaf_queries[leaf].add((name,top_vote))

    print (len(leaf_queries))
    subtree_dict = dict()
    subtree_leaf_label_dict = dict()
    nbr_of_queries = len(q_dict)
    print (len(query_votes_dict), nbr_of_queries)
    most_common_index = 0
    
    # generate the set of subtrees
    while len(query_votes_dict) > 0:
        (seed_label, node_votes) = lf_votes.most_common(most_common_index+1)[most_common_index]
        
        node_y = leaf_dict[seed_label]
        labels = utils.subtree_nodes_with_edge_length(tree, node_y, n)
        subtree = tree.extract_tree_with(labels)
        label_set = set(labels)

        queries_by_subtree = set()
        subtree_query_set = set()

        #gather any other queries that can be used with this subtree
        for label in labels:
            leaf_queries_remove_set = set()
            if label in leaf_queries:
                    
                for leaf_query, top_vote in leaf_queries[label]:
                
                    if leaf_query not in query_votes_dict:
                        leaf_queries_remove_set.add((leaf_query, top_vote))
                        continue
                        
                    if top_vote:
                        subtree_query_set.add(leaf_query)
                        leaf_queries_remove_set.add((leaf_query, top_vote))
                    
                leaf_queries[label].difference_update(leaf_queries_remove_set)
        queries_by_subtree.update(subtree_query_set)

        if len(queries_by_subtree)> 0:
            subtree_dict[subtree] = (seed_label, queries_by_subtree)
            subtree_leaf_label_dict[subtree] = subtree.label_to_node(selection='leaves')
            print ("{} queries in subtree".format(len(queries_by_subtree)))
        # votes_b4 = len(list(lf_votes.elements()))
        
        for query in queries_by_subtree:
            if query in query_votes_dict:
                lf_votes.subtract(query_votes_dict[query])
                query_votes_dict.pop(query)
        #if len(queries_by_subtree)> 0:        
        #    print ("votes before: {}, votes after: {}".format(votes_b4, len(list(lf_votes.elements()))))
        #    print ("queries left: {}".format(len(query_votes_dict)))
        if len(queries_by_subtree) == 0:
            most_common_index += 1
        else:
            most_common_index = 0;
            
    
    # re-assign subtrees to queries
    new_subtree_dict = dict()
    for query, closest_label in query_top_vote_dict.items():
   
        best_subtree = None
        best_distance = 99999999999999999
        for subtree, value in subtree_dict.items():
            leaf_label_dict = subtree_leaf_label_dict[subtree]
            seed_label, _ = value
            if closest_label in leaf_label_dict:
                distance = subtree.distance_between(leaf_label_dict[closest_label], leaf_label_dict[seed_label])
                if distance < best_distance:
                   best_subtree = subtree
                   best_distance = distance
        if best_subtree in new_subtree_dict:
            new_subtree_dict[best_subtree].append(query)
        else:
            new_subtree_dict[best_subtree] = [query]
        
            

    print ('{} seconds re-assigning subtrees'.format(time.perf_counter() - t0))
    final_subtree_count = 0
    
    # ensure there are only max_query_size queries running on a subtree at a time!
    subtree_dict = dict()
    for subtree, query_list in new_subtree_dict.items():
        index = 0
        while index*max_query_size <= len(query_list):
            s_ind = index*max_query_size
            e_ind = (index+1)*max_query_size
            tmp_query_list = query_list[s_ind:e_ind]
            subtree_dict[(subtree, index)] = tmp_query_list
            index += 1
        s_ind = index*max_query_size
        tmp_query_list = query_list[s_ind:]
        subtree_dict[(subtree, index)] = tmp_query_list
    

    # process each subtree query set to get a subalignment
    subalignment_paths = list()
    for subtree_index, query_list in subtree_dict.items():
        subtree, index = subtree_index 

        if len(query_list) == 0:
            continue

        final_subtree_count += 1
        
        tmp_aln_path = "tmp{}/aln".format(run) + ".fa"
        tmp_qaln_path = "tmp{}/qaln".format(run) + "q.fa"
        tmp_output_path = "tmp{}/".format(run) + "/mafft_add_result_{}".format(final_subtree_count)
        tmp_dir = "tmp{}/".format(run)
        try:
            os.mkdir(tmp_dir)
        except OSError as error:
            pass

        # write the temp alignments
        tmp_leaf_dict = subtree.label_to_node(selection='leaves')
        if '' in tmp_leaf_dict:
            del tmp_leaf_dict['']
        tmp_ref_dict = {label : ref_dict[label] for label in tmp_leaf_dict.keys()}
        utils.write_fasta(tmp_aln_path, tmp_ref_dict)
        
        tmp_q_dict = {name : q_dict[name] for name in query_list}
        utils.write_fasta(tmp_qaln_path, tmp_q_dict, aligned=False)
        

        # run mafft-linsi
        # TO DO: need to parameterize the location (add configs and log)
        cmd = '/usr/bin/mafft-linsi --localpair --maxiterate 1000 --quiet --thread {} --add {} {} > {}'.format(16, tmp_qaln_path, tmp_aln_path, tmp_output_path)
        print ('[MAFFT-add] Command used: {}'.format(cmd))
        os.system(cmd)
        
        subalignment_paths.append(tmp_output_path)
        
    print ('{} seconds running MAFFT-add'.format(time.perf_counter() - t0))
    
    merger.mergeAlignments(output_path, subalignment_paths, backbone_aln_path)

    print ('{} seconds merging subalignments'.format(time.perf_counter() - t0))
    print ('Final number of subtrees used:', final_subtree_count)
    #shutil.rmtree("tmp{}".format(run))
    


def parseArgs():
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", type=str,
                        help="Path to reference tree with estimated branch lengths", required=True, default=None)
    
    parser.add_argument("-o", "--output_path", type=str,
                        help="Output path file name", required=False, default="BSCAMPP-EMMA")
    
    parser.add_argument("-c", "--temp_aln_path", type=str,
                        help="Path for query and reference sequence alignment in fasta format", required=True, default=None)

    parser.add_argument("-a", "--backbone_aln_path", type=str,
                        help="Path backbone sequence alignment in fasta format", required=True, default=None)         

    parser.add_argument("-b", "--subtreesize", type=int,
                        help="Integer size of the subtree",
                        required=False, default=2000)
 
    parser.add_argument("-n","--tmpfilenbr", type=int,
                        help="tmp file number",
                        required=False, default=0)
 
    parser.add_argument("-f", "--fragmentflag", type=str2bool,
                        help="boolean, True if queries contain fragments",
                        required=False, default=True)
   
    parser.add_argument("-q", "--queries_path", type=str,
                        help="Path to query sequences in fasta format (backbone alignment separate)",
                        required=False, default="")                        
               
    parser.add_argument("-m", "--model", type=str,
                        help="Model used for edge distances",
                        required=False, default="GTR")
    
    parser.add_argument("-M", "--max_query_size", type=int,
                        help="Maximum number of queries for mafft-linsi-add in a subtree",
                        required=False, default=200)
    
    parser.add_argument("-V", "--votes", type=int,
                        help="Integer number of votes per query sequence",
                        required=False, default=5)
    
    
    parser.add_argument("-v", "--version", action="version", version="1.0.0", help="show the version number and exit")
                       
    return parser.parse_args()


def str2bool(b):
    if isinstance(b, bool):
       return b
    if b.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif b.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')        

if __name__ == "__main__":
    main(parseArgs())
