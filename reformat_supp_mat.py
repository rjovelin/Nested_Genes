# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 12:16:34 2016

@author: RJovelin
"""

# use this function to reformat the supp_table text files into a excel file
def reformat_tables(filename, outputfile):
    '''
    (file) -> file
    Take the supp_table.txt file with list of variants affecting host and
    nested genes and save the file content into the outputfile with a single
    variant per line
    '''
    
    # open file for reading
    infile = open(filename, 'r')
    # get header
    header = infile.readline().rstrip().split()
    
    # modify header
    # remove variation field
    header.remove(header[-1])
    # add variants_ID, variant_start, variant_end, variant_source fields
    fields = ['variant_ID', 'variant_start', 'variant_end', 'variant_source']
    header.extend(fields)
    
    # create a dict to store the data
    data = {}
    
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # use the line number as key
            data[int(line[0])] = line[1:]
    # close file after reading
    infile.close()
    
    # open file for writing
    newfile = open(outputfile, 'w')
      
    # write header to file
    newfile.write('\t'.join(header) + '\n')
    
    # write data to file
    # create a list of keys
    nums = [i for i in data]
    nums.sort()
    # create a counting variable
    j = 1
    # loop over line number in list
    for i in nums:
        # create a variable with variant info
        var = data[i][-1]
        # remove variant info from list
        data[i].remove(data[i][-1])
        # check if multiple variants for the given nested gene
        if ';' in var:
            # multiple variants are recorded
            # get info for each variant
            var = var.split(';')
            # loop over variants
            for k in var:
                # write info about host and nested genes
                newfile.write(str(j) + '\t' + '\t'.join(data[i]) + '\t')
                # write variant info
                newfile.write('\t'.join(k.split(':')) + '\n')
                # update j
                j += 1
        else:
            # a single variant is present
            # parse the variant info, and write the single line
            var = var.split(':')
            data[i].extend(var)
            newfile.write(str(j) + '\t' + '\t'.join(data[i]) + '\n')
            # update j
            j += 1
            
    # close file after writing
    newfile.close()
            
        
    
    
    