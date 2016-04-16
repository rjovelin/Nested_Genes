# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 11:44:35 2015

@author: Richard
"""

import pickle
from elegans_nested_genes import mutations
from elegans_nested_genes import coding_sequences
from elegans_nested_genes import CDS_host_genes
from elegans_nested_genes import variants_host_coding


def count_host_nested_genes(host_genes_pickled, linkage_pickled):
    '''
    (file, file) -> int, dict
    Unpickle the pickled host_genes file containing host: nested gene pairs
    and unpickled the linkage file containing coordinates about all the genes
    in the genome. Return the count of host genes and a dictionnary of nested
    gene counts broken down by nested gene biotype
    '''
    
    # open file for loading
    host_file = open(host_genes_pickled, 'rb')
    # unpickle open file
    host_genes = pickle.load(host_file)
    # close file
    host_file.close()
    
    # open file for loading
    linkage_file = open(linkage_pickled, 'rb')
    # unpickle open file
    linkage = pickle.load(linkage_file)
    # close file
    linkage_file.close()
    
    # get the number of host genes
    N_host_genes = len(host_genes)
    
    # make a set of nested genes
    nested = set()
    # add the nested genes to set
    for gene in host_genes:
        # add all the nested genes
        for nested_gene in host_genes[gene]:
            nested.add(nested_gene)
    
    # create a dictionnary to store the nested gene biotype counts
    nested_biotype = {}
    
    # loop over nested genes and count biotype
    for nested_gene in nested:
        # loop over chromosomes, find the nested gene
        for chromo in linkage:
            for gene in linkage[chromo]:
                if gene == nested_gene:
                    # get biotype
                    biotype = linkage[chromo][gene][-1]
                    # update dictionnary
                    if biotype in nested_biotype:
                        nested_biotype[biotype] += 1
                    else:
                        nested_biotype[biotype] = 1
                    # no need to continue searching if nested_gene is found
                    break
                
    return N_host_genes, nested_biotype


def proportion_of_nesting_in_functional_classes(host_genes_pickled, linkage_pickled):
    '''
    (file, file) -> dict
    Unpickle the files containing information about nested genes and
    coordinates of all genes and return a dictionnary with the proportion
    of nesting for each functional class
    '''
    
    # get the counts of each functional class
    all_gene_file = open(linkage_pickled, 'rb')
    linkage = pickle.load(all_gene_file)
    all_gene_file.close()
    
    # create a dictionnary to store the count of each functional class
    functional_class = {}
    for chromo in linkage:
        for gene in linkage[chromo]:
            biotype = linkage[chromo][gene][-1]
            if biotype in functional_class:
                functional_class[biotype] += 1
            else:
                functional_class[biotype] = 1
                
    # get the functional class count for al nested genes
    nested_genes = count_host_nested_genes(host_genes_pickled, linkage_pickled)[1]
    
    # create dictionnary to store the proportion of nesting for each gene class
    proportions = {}
    for function in nested_genes:
        proportions[function] = round((nested_genes[function] / functional_class[function]) * 100, 2)
    
    return proportions



def proportion_of_nested_genes_by_chromosome(host_genes_pickled, linkage_pickled):
    '''
    (file, file) -> dict
    Unpickle the files containing information about nested genes and
    coordinates of all genes and return a dictionnary with the proportion
    of nesting for each chromsome
    '''
    
    # unpickle files
    all_gene_file = open(linkage_pickled, 'rb')
    linkage = pickle.load(all_gene_file)
    all_gene_file.close()
    
    host_gene_file = open(host_genes_pickled, 'rb')
    host_genes = pickle.load(host_gene_file)
    host_gene_file.close()
    
    # make a set of nested genes
    nested = set()
    # add the nested genes to set
    for gene in host_genes:
        # add all the nested genes
        for nested_gene in host_genes[gene]:
            nested.add(nested_gene)
    
    # create a dictionnary to store counts on each chromosome
    chromosomes = {}
    # loop over all nested genes
    for nested_gene in nested:
        # find the chromosome of the nested gene
        for chromo in linkage:
            for gene in linkage[chromo]:
                if nested_gene == gene:
                    # populate dict
                    if chromo in chromosomes:
                        chromosomes[chromo] += 1
                    else:
                        chromosomes[chromo] = 1
    return chromosomes
    
   
               
def host_nested_variant_relationship(host_nested_variants_file):
    '''
    (file) -> dict
    Read the file containing the variants which positions overlap with the CDS
    of the host gene and the also the nested gene and return a dictionnary 
    of host genes as key and dictionnaries of nested_genes : variants list pairs
    '''
    
    # create dictionnary
    # {host_gene: {nested_gene: [variant1, variant2]}}
    host_variants = {}
    
    # open file for reading
    infile = open(host_nested_variants_file, 'r')
    # read header
    header = infile.readline()
    # go through file and populate dict
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get host gene
            host_gene = line[2]
            # get nested gene
            nested_gene = line[7]
            # get all variants
            variants = line[-1]
            # parse variants to get the variant IDs
            variants = variants.split(';')
            for i in range(len(variants)):
                variants[i] = variants[i].split(':')
            # create list to get the variants
            variant_IDs = []
            for i in range(len(variants)):
                for item in variants[i]:
                    if item.startswith('WBVar'):
                        variant_IDs.append(item)
            # populate dict
            if host_gene in host_variants:
                # check if nested_gene is key of inner dict
                if nested_gene in host_variants[host_gene]:
                    host_variants[host_gene][nested_gene].extend(variant_IDs)
                else:
                    host_variants[host_gene][nested_gene] = variant_IDs
            else:
                host_variants[host_gene] = {}
                host_variants[host_gene][nested_gene] = variant_IDs
    
    # close file after reading
    infile.close()
    
    return host_variants


def get_nested_genes(host_nested_variants_file):
    '''
    (file) -> set
    Read the file containing the variants which positions overlap with the CDS
    of the host gene and the also the nested gene and return a set of nested genes    
    '''
    
    # create a dict of host genes keys and dicts of nested gene : variant list pairs
    host_variants = host_nested_variant_relationship(host_nested_variants_file)
    
    # create a set of nested genes
    nested_genes = set()
    for gene in host_variants:
        for nested in host_variants[gene]:
            nested_genes.add(nested)
    
    return nested_genes
  
    
def get_nested_variants(host_nested_variants_file):
    '''
    (file) -> set
    Read the file containing the variants which positions overlap with the CDS
    of the host gene and the also the nested gene and return the set of
    variant IDs
    '''
    
    # create a dict of host genes keys and dicts of nested gene : variant list pairs
    host_variants = host_nested_variant_relationship(host_nested_variants_file)
    
    # create a set of variants
    nested_variants = set()
    for gene in host_variants:
        for nested in host_variants[gene]:
            for var in host_variants[gene][nested]:
                nested_variants.add(var)
    return nested_variants    

    
def count_all_variants(host_nested_variants_file):
    '''
    (file) -> int, int, int
    Read the file containing the variants which positions overlap with the CDS
    of the host gene and the also the nested gene and return the count of 
    affected host genes, affected nested genes and count of variants
    '''
    
    # create a dict of host genes keys and dicts of nested gene : variant list pairs
    host_variants = host_nested_variant_relationship(host_nested_variants_file)
    
    # get the number of host genes
    N_host = len(host_variants)
    
    # get the number of nested genes:
    # create a set of nested genes
    nested_genes = get_nested_genes(host_nested_variants_file)    
    N_nested = len(nested_genes)

    # get the number of variants
    nested_variants = get_nested_variants(host_nested_variants_file)
    N_variants = len(nested_variants)
    
    return N_host, N_nested, N_variants
    
    
    
def variants_coord(gff_file):
    '''
    (file) -> dict
    Return a dictionnary of variants: coordinates pairs for all variants
    in the gff file. Coordinates is a list including the chromosome, start and
    end positions, source and mutation type
    '''
    
    # create a dict to store the variants : coord pairs
    # {var_ID: [chromo, start, end, source, type]}
    variants = {}
        
    # open file for reading
    infile = open(gff_file, 'r')
        
    # loop over file
    for line in infile:
        line = line.rstrip()
        # check that variant in line
        if 'WBVar' in line:
            # make a list
            line = line.split()
            # get chromosome
            chromo = line[0]
            # get variation source
            source = line[1]
            # variant name
            var_name = line[8][line[8].index('WBVar'): line[8].index(';', line[8].index('WBVar'))]
            # get start and end positions
            start, end = int(line[3]), int(line[4])
            # mutation type
            var_type = line[2]
            variants[var_name] = [chromo, start, end, source, var_type]

    # close file
    infile.close()
    
    return variants
    
           
def cloned_genes(gene_ID_file):
    '''
    (file) -> dict
    Return a dictionnary of gene ID: gene name for genes that have been cloned
    '''

    # create a dict of gene ID: gene name pairs
    genes = {}
    
    # open file for reading
    infile = open(gene_ID_file, 'r')
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split(',')
            # use only live genes
            if line[-1] == 'Live':
                # get the gene ID
                geneID = line[1]
                # cloned name has '-'
                if '-' in line[2]:
                    gene_name = line[2]
                    genes[geneID] = gene_name
    # close file
    infile.close()
    
    return genes
    
   
def host_contained_variants_to_file(host_nested_variants_file, outputfile):
    '''
    (file, file) -> file
    Read the file containing the variants which positions overlap with the CDS
    of the host gene and also the nested gene and copy content to the 
    outputfile with variants that are entirely within the host gene
    '''
    
    # open file for writing
    newfile = open(outputfile, 'w')
    
    # open file for reading
    infile = open(host_nested_variants_file, 'r')
    # read header
    header = infile.readline()
    
    # write header to newfile
    newfile.write(header)  
    
    # set new counter
    j = 1
        
    # go through infile
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get the coordinates of the host gene
            hoststart, hostend = int(line[3]), int(line[4])
            # parse the variants
            line[-1] = line[-1].split(';')
            # parse each variant to get the coordinates
            for i in range(len(line[-1])):
                line[-1][i] = line[-1][i].split(':')
            # create a list of variant to remove
            to_remove = []
            # add variant to to_remove if variant coord are not within host coord
            for i in range(len(line[-1])):
                varstart, varend = int(line[-1][i][1]), int(line[-1][i][2])
                if varstart < hoststart or varend > hostend:
                    to_remove.append(line[-1][i])
            for i in range(len(to_remove)):
                line[-1].remove(to_remove[i])
            # copy line to newfile if variants remain
            if len(line[-1]) != 0:
                # change counter
                line[0] = str(j)
                # convert back list to strings
                for i in range(len(line[-1])):
                    line[-1][i] = ':'.join(line[-1][i])
                line[-1] = ';'.join(line[-1])
                # copy line to newfile
                for item in line[:-1]:
                    newfile.write(item + '\t')
                newfile.write(line[-1] + '\n')
                # update counter
                j += 1
    
    # close files
    infile.close()
    newfile.close()



def non_overlapping_variants(host_nested_variants_file, linkage_pickled, gff_file):
    '''
    (file, file, file) -> dict   
    Return a dictionnary with host genes as key and a list of variants affecting
    the host's CDS but not the sequence of nested gene that are contained within
    the host gene
    '''
    
    # get the genes' coordinates
    # {chromo: {gene:[chromosome, start, end, sense, biotype]}}
    linkagefile = open(linkage_pickled, 'rb')
    linkage = pickle.load(linkagefile)
    linkagefile.close()
    
    # create a dict with host genes as keys
    # {host_gene: {nested_gene: [variant1, variant2]}}
    host_genes = host_nested_variant_relationship(host_nested_variants_file)
       
    # create a set of variants affecting the nested gene and the CDS of the host
    host_nested_variants = set()
    for gene in host_genes:
        for nested_gene in host_genes[gene]:
            for var in host_genes[gene][nested_gene]:
                host_nested_variants.add(var)
       
    # create a dict with positions of indel variants
    # {LG: {WBVAr: [start, end, variation_type]}}
    variants = mutations(gff_file)
    
    # create a dict with the CDS positions
    # {chromo: {gene:[(CDS_start, CDS_end), (CDS_start, CDS_end)]}}
    coding = coding_sequences(gff_file)
    
    # create a dict with the CDS positions of the host_genes
    # {chromo: {gene:[(CDS_start, CDS_end), (CDS_start, CDS_end)]}}
    host_genes_coding = CDS_host_genes(coding, host_genes)
    
    # create a dictionnary of variants affecting the CDS of host but not the nested genes
    # { host_gene: [variant1, variant2]}    
    host_variants_only = {}    
      
    # loop over host_genes_coding
    for chromo in host_genes_coding:
        # loop over host genes on chromo
        for gene in host_genes_coding[chromo]:
            # loop over variants on chromo
            for var in variants[chromo]:
                # check that var does not affect host and nested genes
                if var not in host_nested_variants:
                    # check that var is contained within the host gene
                    host_start = linkage[chromo][gene][1]
                    host_end = linkage[chromo][gene][2]
                    var_start = variants[chromo][var][0]
                    var_end = variants[chromo][var][1]
                    # check variant starting position is within the host gene
                    if var_start >= host_start and var_start <= host_end:
                        # check that variant is contained witin the host gene
                        host_coord = set(range(host_start, host_end +1))
                        var_coord = set(range(var_start, var_end + 1))
                        if var_coord.issubset(host_coord):
                            # check that variant affect the host CDS
                            for CDS in host_genes_coding[chromo][gene]:
                                # get the coordinates of each CDS
                                CDS_coord  = set(range(CDS[0], CDS[1] + 1))
                                # check if the gene's CDS overlap with variant coordinates
                                if len(var_coord.intersection(CDS_coord)) != 0:
                                    # populate dict
                                    if gene in host_variants_only:
                                        host_variants_only[gene].append(var)
                                    else:
                                        host_variants_only[gene] = [var]
                                     # no need to evaluate all CDS if CDS is affected by variant
                                    break
                                
    return host_variants_only                    
                    

def alternative_variants(host_nested_variants_file, linkage_pickled, gff_file, outputfile):
    '''
    (file, file, file, file) -> file
    Save to file the list of variants affecting the coding sequence of host genes
    but not the sequence of the nested genes
    '''
    
    # create a dictionnary of variants affecting the CDS of host but not the nested genes
    # { host_gene: [variant1, variant2]}    
    host_variants_only = non_overlapping_variants(host_nested_variants_file, linkage_pickled, gff_file)
    
    # create a dict to store the variants : coord pairs
    # {var_ID: [chromo, start, end, source, type]}
    variants = variants_coord(gff_file)
    
    # unpickle dict linkage
    linkagefile = open(linkage_pickled, 'rb')
    linkage = pickle.load(linkagefile)
    linkagefile.close()
    
    # open file for writing
    newfile = open(outputfile, 'w')
    # write header
    header = '\t'.join(['Variant', 'Chromosome', 'Variant_start', 'Variant_end', 'Host_gene', 'Host_start', 'Host_end', 'Host_sense'])
    newfile.write(header + '\n')
    
    # loop over dict and save content to file
    for gene in host_variants_only:
        # find the coordinates of the host_gene
        for LG in linkage:
            for locus in linkage[LG]:
                if locus == gene:
                    chromo = LG
                    host_start = linkage[LG][gene][1]
                    host_end = linkage[LG][gene][2]
                    sense = linkage[LG][gene][3]
                    break
        # get the coordinates of each variant   
        for variant in host_variants_only[gene]:
            var_start = variants[variant][1] 
            var_end = variants[variant][2]
            # write information to file
            info = '\t'.join([variant, chromo, str(var_start), str(var_end), gene, str(host_start), str(host_end), sense]) 
            newfile.write(info + '\n')
    
    # close file after writing
    newfile.close()

   
    
def nested_variant_phenotypes(nested_alleles_phenotypes_file):
    '''
    (file) -> dict
    Take the file containing the phenotypes for the alleles of the host genes
    and return a dictionnary of host_genes as keys and dictionnaries of variants:
    list of phenotypes pairs as values
    '''
    
    # open file for reading
    infile = open(nested_alleles_phenotypes_file, 'r')
    # skip header
    infile.readline()
    
    # create a dictionnary
    # {host: variant: {pheno1, pheno2}}
    pheno = {}
    
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # get the host gene ID
            host = line[0]
            # get the variant ID
            variant = line[3]
            # get the description of the phenotype
            phenotype = line[5]
            # populate the dict
            # check that host gene is key to dict
            if host in pheno:
                # check that variant is key to inner dict
                if variant in pheno[host]:
                    pheno[host][variant].add(phenotype)
                else:
                    pheno[host][variant] = set(phenotype)
            else:
                pheno[host] = {}
                pheno[host][variant] = set(phenotype)
                
    # close file ater reading
    infile.close()
    
    return pheno

def alleles_phenotypes(all_genes_phenotypes_file, linkage_pickled):
    '''
    (file, file) -> dict
    Return a dictionnary with protein coding gene as key and an inner
    dictionnary as value with variants: list of phenotypes pairs
    Restrict the genes to protein coding genes using the linkage file
    '''
    
    # create a dictionnary of {chromo: {gene: [coordinates]}}
    # open file for unpickling
    linkage_file = open(linkage_pickled, 'rb')
    # load dictionnary from file
    linkage = pickle.load(linkage_file)
    # close file
    linkage_file.close()
    
    # open file for reading
    infile = open(all_genes_phenotypes_file, 'r')
    # skip header
    infile.readline()
    
    # create dict {gene: var: [pheno1, pheno2]}
    alleles = {}
    
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            if '"' in line:
                line = line.replace('"', '')
            if ',' in line:
                line = line.split(',')
            else:
                line = line.split()
            for item in line:
                if item.startswith('WBGene'):
                    gene = item
                elif item.startswith('WBVar'):
                    var = item
                elif item.startswith('WBPheno'):
                    pheno = item
            # populate dict
            if gene in alleles:
                if var in alleles[gene]:
                    alleles[gene][var].append(pheno)
                else:
                    alleles[gene][var] = [pheno]
            else:
                alleles[gene] = {}
                alleles[gene][var] = [pheno]
            
    #close file after reading
    infile.close()
    
    # remove non protein coding genes
    to_remove = []
    for gene in alleles:
        for chromo in linkage:
            if gene in linkage[chromo]:
                if linkage[chromo][gene][-1] != 'protein_coding':
                    to_remove.append(gene)
    for gene in to_remove:
        del alleles[gene]
        
    return alleles

    
    
    
def host_with_phenotype_causing_vartiants(all_genes_phenotypes_file, linkage_pickled, host_nested_variants_file):
    '''
    (file, file, file) -> dict
    From the file with variants contained entirely within the host protein 
    coding genes and affecting the coding sequence of the host and the sequence
    of nested genes, filter out genes that do not have a phenotype using the
    all_genes_phenotypes and linkage_pickled files and return a dictionnary of
    host protein coding genes with dicts of nested_genes and list of
    phenotype-causing variants as value
    '''
    
    # create a dict {gene: var: [pheno1, pheno2]} for all genes with phenotype
    alleles = alleles_phenotypes(all_genes_phenotypes_file, linkage_pickled)
    
    # create a dict {host_gene: {nested_gene: [variant1, variant2]}} 
    # for host genes with variants and nested genes
    hosts = host_nested_variant_relationship(host_nested_variants_file)
    
    # create a dict {host_gene: {nested_gene: [variant1, variant2]}}
    # for hosts with phenotype-causing variants also affecting the nested gene
    pheno = {}
    for gene in hosts:
        if gene in alleles:
            pheno[gene] = hosts[gene]
    
    # remove variants that do not cause a phenotype
    # create a set phenotype-causing variants
    variants = set()
    for gene in alleles:
        for var in alleles[gene]:
            variants.add(var)
    
    # loop over host genes, remove non-phenotype causing alleles
    for gene in pheno:
        for nested_gene in pheno[gene]:
            # create a list of variants to remove:
            to_remove = []
            # add var to list if no associated phenotype
            for var in pheno[gene][nested_gene]:
                if var not in variants:
                    to_remove.append(var)
            # remove non-phenotype causing variants
            for var in to_remove:
                pheno[gene][nested_gene].remove(var)
                
    # create a dict with host genes as keys and dicts of nested_gene : list
    # of variants affecting the host's CDS and the sequence of nested gene
    # and causing a phenotpe
    pheno_var = {}
    
    # loop over genes in pheno
    for gene in pheno:
        # loop over nested_genes
        for nested_gene in pheno[gene]:
            # keep nested_gene if it is affected by some variants
            if len(pheno[gene][nested_gene]) != 0:
                if gene in pheno_var:
                    pheno_var[gene][nested_gene] = list(pheno[gene][nested_gene])
                else:
                    pheno_var[gene] = {}
                    pheno_var[gene][nested_gene] = list(pheno[gene][nested_gene])
    
    return pheno_var




def host_with_phenotype_to_file(all_genes_phenotypes_file, linkage_pickled, host_nested_variants_file, outputfile):
    '''
    (file, file, file, file) -> file
    From the file with variants contained entirely within the host protein 
    coding genes and affecting the coding sequence of the host and the sequence
    of nested genes, select the host protein coding genes that have a phenotype
    using the all_genes_phenotype_file and linkage_picked files and write to
    outputfile information about the host genes, the phenotype-causing variants
    and the affected nested genes
    '''
    
    # create a dict of hosts and nested genes affected by phenotype-causing variants
    pheno_var = host_with_phenotype_causing_vartiants(all_genes_phenotypes_file, linkage_pickled, host_nested_variants_file)
    # create a set of nested genes affected by the phenotype-causing variants
    altered_nested = set()
    for gene in pheno_var:
        for nested_gene in pheno_var[gene]:
            altered_nested.add(nested_gene)
    # create a set of phenotype-causing variants
    variants = set()
    for gene in pheno_var:
        for nested_gene in pheno_var[gene]:
            for var in pheno_var[gene][nested_gene]:
                variants.add(var)
          
    # open file for writing
    newfile = open(outputfile, 'w')
    
    # write header
    header = '\t'.join(['#', 'chromosome', 'host_gene', 'host_start',
                        'host_end', 'host_sense', 'host_biotype', 'nested_gene',
                        'nested_start', 'nested_end', 'nested_sense',
                        'nested_biotype', 'variation'])
    newfile.write(header + '\n')
    
    # set new counter
    j = 1
    
    # open file with host, nested genes and variant information
    infile = open(host_nested_variants_file, 'r')
    # skip header
    infile.readline()
    
    # go through infile
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get host gene
            host_gene = line[2]
            nested_gene = line[7]
            # check that host_gene is in pheno_var and that its nested_gene
            # is affected by the variant
            if host_gene in pheno_var and nested_gene in altered_nested:
                # parse the variants
                line[-1] = line[-1].split(';')
                # parse each variant to get variant ID
                for i in range(len(line[-1])):
                    line[-1][i] = line[-1][i].split(':')
                # create a list of variant to remove
                to_remove = []
                # add variant to to_remove if variant_ID not causing phenotype
                for i in range(len(line[-1])):
                    var_ID = line[-1][i][0]
                    if var_ID not in variants:
                        to_remove.append(line[-1][i])
                # remove variants not causing a phenotype
                for i in range(len(to_remove)):
                    line[-1].remove(to_remove[i])
                # copy to newfile if variants for the nested gene
                if len(line[-1]) != 0:
                    # change counter
                    line[0] = str(j)
                    # convert back list to strings
                    for i in range(len(line[-1])):
                        line[-1][i] = ':'.join(line[-1][i])
                    line[-1] = ';'.join(line[-1])
                    # copy line to newfile
                    for item in line[:-1]:
                        newfile.write(item + '\t')
                    newfile.write(line[-1] + '\n')
                    # update counter
                    j += 1
    
    # close files
    infile.close()
    newfile.close()

        
    
def fraction_nested_genes_altered(host_nested_variants_file):
    '''
    (file) -> dict
    Return a dictionnary with variant ID as key and a list containing the 
    the proportions of the nested gene length altered by the variant
    '''
    
    # open file for reading
    infile = open(host_nested_variants_file, 'r')
    # skip header
    infile.readline()
    
    # create dict {variant_ID: [num]}
    nested_alteration = {}
        
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get variant ID
            nested_gene = line[7]            
            # get the nested gene length
            nested_start, nested_end = int(line[8]), int(line[9])
            nested_length = nested_end - nested_start + 1
            nested_coord = set(range(nested_start, nested_end + 1))
            # parse the variant information
            # parse the multiple variant
            line[-1] = line[-1].split(';')
            # get the coorrdinates the each variant
            for i in range(len(line[-1])):
                line[-1][i] = line[-1][i].split(':')
            # get variant ID, variant coordinates for all variants
            for i in range(len(line[-1])):
                # get the variant ID
                variant = line[-1][i][0]
                # get the variant coordinates
                variant_start, variant_end = int(line[-1][i][1]), int(line[-1][i][2])
                variant_coord = set(range(variant_start, variant_end + 1))
                # compute  fraction of nested gene length altered
                fraction = (len(variant_coord.intersection(nested_coord)) / nested_length) * 100
                # populate dict
                if variant in nested_alteration:
                    nested_alteration[variant].append(fraction)
                else:
                    nested_alteration[variant] = [fraction]
    
    # close file after reading
    infile.close()

    return nested_alteration                    
                    
                
            
            
            


    
    
    
    
    
        
    