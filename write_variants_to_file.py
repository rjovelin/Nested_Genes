# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 15:08:46 2015

@author: Richard
"""

def variants_to_file(outputfile, host_nested_variants, linkage):
    '''
    (file) -> file
    Find the variants altering both the host genes CDS and the sequence of the
    nested gene and save information about the host, nested genes and variants
    to an outputfile
    '''
       
    # open file for writing
    newfile = open(outputfile, 'w')
    
    # write header
    header = '\t'.join(['#', 'chromosome', 'host_gene', 'host_start', 'host_end',
                        'host_sense', 'host_biotype', 'nested_gene',
                        'nested_start', 'nested_end', 'nested_sense', 
                        'nested_biotype', 'variation'])
    newfile.write(header + '\n')

    # set counter    
    j = 1
    
    # loop over the host genes in host_nested_variants
    for host_gene in host_nested_variants:
        # get the coordinates of the host gene
        for LG in linkage:
            for gene in linkage[LG]:
                if gene == host_gene:
                    host_chromo = linkage[LG][host_gene][0]
                    host_start = str(linkage[LG][host_gene][1])
                    host_end = str(linkage[LG][host_gene][2])
                    host_sense = linkage[LG][host_gene][3]
                    host_biotype = linkage[LG][host_gene][4]
                    break
        # loop over the nested genes
        for nested_gene in host_nested_variants[host_gene]:
            # get the coordinates of the nested_gene
            for LG in linkage:
                for gene in linkage[LG]:
                    if gene == nested_gene:
                        nested_start = str(linkage[LG][nested_gene][1])
                        nested_end = str(linkage[LG][nested_gene][2])
                        nested_sense = linkage[LG][nested_gene][3]
                        nested_biotype = linkage[LG][nested_gene][4]
                        break
            # loop over all variants affecting the nested gene
            # initialize string to store variation information
            variation = ''
            for variant in host_nested_variants[host_gene][nested_gene]:
                # convert elements of variant to string
                variant = [str(i) for i in variant]
                # join all features of a single variant by ':'
                variant = ':'.join(variant)
                variation += (variant + ';')
            # remove last '_'
            variation = variation[:-1]
            
            # write host gene and combined variation for each nested gene to file
            nested_info = '\t'.join([str(j), host_chromo, host_gene, host_start, host_end,
                           host_sense, host_biotype, nested_gene, nested_start,
                           nested_end, nested_sense, nested_biotype, variation])
            newfile.write(nested_info + '\n')
            
            # update j
            j += 1
            
    # close file after writing
    newfile.close()