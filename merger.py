# this is adapted from EMMA at https://github.com/c5shen/EMMA/
import time, os
from alignmentTools import Alignment, CompactAlignment, \
        compact, ExtendedAlignment

def mergeAlignments(output_path, subalignment_paths, backbone_path):
    #Configs.log('Merging all sub-alignments with transitivity with ' \
    #        'singletons from queries collapsed')
    #start = time.time()

    #aln_dir = os.path.join(Configs.outdir, 'sub-alignments')

    init_aln = Alignment(); init_aln.read_file_object(backbone_path)
    new_aln = compact(init_aln)
    for i in range(0, len(subalignment_paths), 1):
        path = subalignment_paths[i]
        subaln = Alignment(); subaln.read_file_object(path)
        new_aln.merge_in(compact(subaln))
        del subaln

    new_aln.write(output_path, 'FASTA')
    #Configs.log('Finished merging all sub-alignments, output file: {}'.format(
    #    output_path))
    #time_merge = time.time() - start
    #Configs.runtime('Time to merge all sub-alignments (s): {}'.format(time_merge))

