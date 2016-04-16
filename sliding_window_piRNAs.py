# -*- coding: utf-8 -*-
"""
Created on Sat May 23 01:59:53 2015

@author: Richard
"""

from mir8205_analysis import convert_elegans_genome
from elegans_nested_genes import *


def sliding_window_piRNA_counts(genome_fasta, linkage_pickled, host_genes_pickled, window_size, outputfile):
    '''
    (file, file, file, int, str) -> dict
    Save to outputfile a table with window number and the counts in each window
    of length window_size of piRNAs located on chromosome IV that are nested
    and not nested length
    '''
    
    # convert genome to dict
    genome = convert_elegans_genome(genome_fasta)
    
    # get the nested genes : host pairs
    nested_genes = nested_host_pairs(host_genes_pickled)
    
    # get the genes' coordinates
    linkage = unpickle_dict(linkage_pickled)
    
    # make lists of size window that contains only 0s:
    # each value in the list is the count of position for the range [0 - window[ etc
    range_counts_nested = [0] * (len(genome['IV']) // window_size)
    range_counts_not_nested = [0] * (len(genome['IV']) // window_size)

    # loop over chromosome IV
    for gene in linkage['IV']:
        # check if gene is piRNA
        if linkage['IV'][gene][-1] == 'piRNA':
            # get starting position
            start = linkage['IV'][gene][1]
            # determine the index in the list range_count where the position should be added
            which_range = start // window_size
            if which_range == len(range_counts_nested):
                which_range -= 1
            # check if piRNA or not
            if gene in nested_genes:
                # update count at this index
                range_counts_nested[which_range] += 1
            else:
                # update count at this index
                range_counts_not_nested[which_range] += 1
    
    # open file for writing
    newfile = open(outputfile, 'w')
    # write header
    newfile.write('Range' + '\t' + 'Lower_point' + '\t' + 'Midpoint' + '\t' + 'Higher_point' + '\t' + 'Count_nested' + '\t' + 'Coun_not_nested' + '\n')
        
    # loop over indices of list
    for i in range(len(range_counts_nested)):
        newfile.write('[' + str(i * window_size) + '-' + str((i * window_size) + window_size -1) + ']' + '\t')
        newfile.write(str(i * window_size) + '\t')
        newfile.write(str(int(((i * window_size) + window_size) / 2)) + '\t')
        newfile.write(str((i * window_size) + window_size) + '\t')
        newfile.write(str(range_counts_nested[i]) + '\t' + str(range_counts_not_nested[i]) + '\n')
        
    newfile.close()

    
    