# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:27:04 2015

@author: Richard
"""

# import functions from modules
from elegans_nested_genes import *
from variants_nested_genes import *
import numpy as np
import math
import pickle



def phenotype_causing_variants_ID(host_nested_variants_file):
    '''
    (file) -> set
    Return a set with the ID of the phenotype-causing mutations that alter
    the host's CDS and the sequence of the nested gene
    '''
    
     # create a dict from the file of host, nested genes and phenotype-causing variants
    # {host_gene: {nested_gene: [variant1, variant2]}}
    hosts = host_nested_variant_relationship(host_nested_variants_file)
    
    # create a set of phenotype-causing variants atering host and nested genes
    phenotype_variants = set()
    for gene in hosts:
        for nested_gene in hosts[gene]:
            for var in hosts[gene][nested_gene]:
                phenotype_variants.add(var)
    
    return phenotype_variants



def variant_mutation_type(gff_file, host_nested_variants_file):
    '''
    (file, file) -> dict
    Return a dictionnary with the mutation type causing a phenotype in the 
    host gene and altering the sequence of nested genes as key and a set of
    corresponding variants as value 
    '''
    
    # create a dict with the variant coordinates and mutation type
    # {var_ID: [chromo, start, end, source, type]}
    variants = variants_coord(gff_file)
    
    # create a set of phenotype-causing variants atering host and nested genes
    phenotype_variants = phenotype_causing_variants_ID(host_nested_variants_file)
    
    # create a dict of mutation_type : {variant1, variant2}
    mutation_type = {}
    for var in phenotype_variants:
        mutation = variants[var][-1]
        if mutation in mutation_type:
            mutation_type[mutation].add(var)
        else:
            mutation_type[mutation] = set()
            mutation_type[mutation].add(var)
    
    return mutation_type


def variant_mutation_consequence(gff_file, host_nested_variants_file):
    '''
    (file, file) -> dict
    Return a dictionnary with the annotated mutation consequence from the gff file
    as key and a set of phenotype-causing mutations that alter the host's CDS
    and the sequence of the nested gene as value
    '''
    
    # create a set of phenotype-causing variants atering host and nested genes
    phenotype_variants = phenotype_causing_variants_ID(host_nested_variants_file)
    
    # open file for reading
    infile = open(gff_file, 'r')
    
    # create a dict to store the effect of mutation
    effect = {}
    
    # loop over file, populate dict
    for line in infile:
        # check that line contains variant information
        if 'WBVar' in line:
            # make a list from line
            line = line.rstrip().split()
            # get the variant ID
            var_name = line[8][line[8].index('WBVar'): line[8].index(';', line[8].index('WBVar'))]
            # check that the variant ID causes a phenotype for the host gene
            if var_name in phenotype_variants:
                # get the mutation effect
                consequence = line[8][line[8].index('consequence=')+12: ]
                # add variant_ID
                if consequence in effect:
                    effect[consequence].add(var_name)
                else:
                    effect[consequence] = set()
                    effect[consequence].add(var_name)
                    
    
    infile.close()                
    return effect

    

def variant_mutation_type_count(gff_file, host_nested_variants_file):
    '''
    (file, file) -> dict
    Return a dictionnary with the mutation type causing a phenotype in the 
    host gene and altering the sequence of nested genes as key and its count as
    value 
    '''
    
    # create a dict with the mutation_type : {variant1, variant2}
    mutation_type = variant_mutation_type(gff_file, host_nested_variants_file)
    
    # create a dict with mutation_type : count pairs
    mutation_type_counts = {}
    
    # populate dict
    for mutation in mutation_type:
        mutation_type_counts[mutation] = len(mutation_type[mutation])
    
    return mutation_type_counts
    
        
    
def variant_mutation_consequence_count(gff_file, host_nested_variants_file):
    '''
    (file, file) -> dict
    Return a dictionnary with the annotated mutation consequence from the gff file
    as key and its count as value for the phenotype-causing mutations that 
    alter the host's CDS and the sequence of the nested gene
    '''
    
    # create a dict with the mutation_effect : {variant1, variant2}
    mutation_effect = variant_mutation_consequence(gff_file, host_nested_variants_file)
    
    # create a dict with mutation_effect : count pairs
    mutation_effect_counts = {}
    
    # populate dict
    for effect in mutation_effect:
        mutation_effect_counts[effect] = len(mutation_effect[effect])
    
    return mutation_effect_counts
    


def function_nested_genes_phenotype_causing_variants(host_nested_variants_file):
    '''
    (file) -> dict
    Return a dictionnary with functional class : count pairs for the nested genes
    altered by phenotype_causing mutations using the file with information about
    host, nedted genes and phenotype-variants
    '''
    
    # open file for reading
    infile = open(host_nested_variants_file, 'r')
    # skip header
    infile.readline()
    
    # create a dict with functional class a key and count as value
    nested_function = {}
    # loop over the file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get the biotype of the nested gene
            biotype = line[11]
            if biotype in nested_function:
                nested_function[biotype] += 1
            else:
                nested_function[biotype] = 1
                
    # close file
    infile.close()
    
    return nested_function
    
    
        
    
def phenotype_causing_variant_lengths(host_nested_variants_file):
    '''
    (file) -> dict
    Return a dictionnary with variant_ID as key and length of variant as value
    '''
    
    # create a dictionnary to store the variant_ID : variant length pairs
    variants = {}
    
    # open file for reading
    infile = open(host_nested_variants_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip().split()
        # parse variants
        line[-1] = line[-1].split(';')
        for i in range(len(line[-1])):
            line[-1][i] = line[-1][i].split(':')
        # go through each variation, get the variant_ID and compute its length
        for i in range(len(line[-1])):
            variant_ID = line[-1][i][0]
            start = int(line[-1][i][1])
            end = int(line[-1][i][2])
            # populate dict
            variants[variant_ID] = end - start + 1
            
    # close file after reading
    infile.close()
    
    return variants
            
        
def mean_phenotype_causing_variant_length(gff_file, host_nested_variants_file, mutation_type_or_effect):
    '''
    (file) -> dict
    Return a dictionnary with the mutation type or the mutation consequence
    as key and a list with the mean variant length and standard error as value
    for phenotype-causing variants
    '''
    
    # select mutation type or mutation effect    
    if mutation_type_or_effect == 'type':
        # get the variant length for the type of mutation
        mutations = variant_mutation_type(gff_file, host_nested_variants_file)
    elif mutation_type_or_effect == 'effect':
        mutations = variant_mutation_consequence(gff_file, host_nested_variants_file)
        
    # create a dict to store the variant ID : variant length pair
    variants = phenotype_causing_variant_lengths(host_nested_variants_file)
    
    # create a dict to store the mutation_type or mutation_effect : variant length pairs
    mutation_length = {}
    
    # loop over the variants
    for var in variants:
        # check which mutation_type or mutation_effect is associated with var
        for key in mutations:
            if var in mutations[key]:
                # update dict with mutation_type or effect and variant length
                if key in mutation_length:
                    mutation_length[key].append(variants[var])
                else:
                    mutation_length[key] = [variants[var]]
    
    # create a dict to store the mutation type or effect as key and the mean
    # variant length as value
    means = {}
    
    for key in mutation_length:
        # compute mean
        means[key] = [np.mean(mutation_length[key])]
        # add standard error of the mean
        means[key].append(np.std(mutation_length[key]) / math.sqrt(len(mutation_length[key])))
        
    return means
        

    
def variants_phenotype_non_host_genes(all_genes_phenotypes_file, linkage_pickled, host_genes_pickled, gff_file):
    '''
    (file, file, file, file) -> dict
    Return a dict with non-host genes as keys as a 2-item list as value 
    containing 2 non-overlapping phenotype-causing variants 
    '''
    # create a dict to store the variants : coord pairs
    # {var_ID: [chromo, start, end, source, type]}
    variants = variants_coord(gff_file)
        
    # create dict with phenotype-causing variants affecting protein coding genes
    # {gene: {var: [pheno1, pheno2]}}
    alleles = alleles_phenotypes(all_genes_phenotypes_file, linkage_pickled)
    
    # get all the genes' coordinates
    linkagefile = open(linkage_pickled, 'rb')
    linkage = pickle.load(linkagefile)
    linkagefile.close()    
    
    # get the host_genes
    hostgenesfile = open(host_genes_pickled, 'rb')    
    host_genes = pickle.load(hostgenesfile)
    hostgenesfile.close()
    
    # create a dict with CDS coordinates for all protein-coding genes
    # {chromo: {gene:[(CDS_start, CDS_end), (CDS_start, CDS_end)]}}  
    coding = coding_sequences(gff_file)
    
    # create a dict with the CDS coordinates of the non-host genes having a phenotype
    # {chromo: {gene:[(CDS_start, CDS_end), (CDS_start, CDS_end)]}}
    non_host_coding = {}
    # loop over each chromosome
    for chromo in coding:
        # loop over all genes on chromo
        for gene in coding[chromo]:
            # check if gene is host gene
            if gene not in host_genes:
                # check that gene has a phenotype
                if gene in alleles:
                    # grab the CDS coordinates
                    if chromo in non_host_coding:
                        non_host_coding[chromo][gene] = list(coding[chromo][gene])
                    else:
                        non_host_coding[chromo] = {}
                        non_host_coding[chromo][gene] = list(coding[chromo][gene])
    
    # create a dictionnary of genes: variants affecting the CDS and haping a phenotype
    # {non_host_gene: [variant1, variant2]}    
    variants_pheno = {}
    
    # loop over chromosomes
    for chromo in non_host_coding:
        # loop over non-host genes having a phenotype      
        for gene in non_host_coding[chromo]:
            # get the gene coordinates
            gene_start = linkage[chromo][gene][1]
            gene_end = linkage[chromo][gene][2]
            gene_coord = set(range(gene_start, gene_end + 1))
            # loop over all variants affecting this gene
            for var in alleles[gene]:
                # check that variant is in the GFF file
                if var in variants:
                    # get the variants coordinates
                    var_start, var_end = variants[var][1], variants[var][2]
                    # check that variant is indel
                    if var_start != var_end:
                        var_coord = set(range(var_start, var_end + 1))
                        # check if the variant is contained within the boundaries of the host gene
                        if var_coord.issubset(gene_coord):
                            # check that the variant affects the CDS
                            for CDS in non_host_coding[chromo][gene]:
                                CDS_start, CDS_end = CDS[0], CDS[1]
                                CDS_coord = set(range(CDS_start, CDS_end + 1))
                                if len(CDS_coord.intersection(var_coord)) != 0:
                                    # populate dict
                                    if gene in variants_pheno:
                                        variants_pheno[gene].append(var)
                                    else:
                                        variants_pheno[gene] = [var]
                                    # exit loop is variant affect at least 1 CDS
                                    break
                        
    # remove genes that have only 1 variant
    to_remove = []
    for gene in variants_pheno:
        if len(variants_pheno[gene]) < 2:
            to_remove.append(gene)
    for gene in to_remove:
        del variants_pheno[gene]
    
    # for each gene, find a pair of non-overlapping variants
    # create dict with non-host gene as key and 2-item list as value
    non_host_genes = {}
    # loop over genes in variants_pheno
    for gene in variants_pheno:
        # set boolean
        not_found = True        
        # loop over variants, compare coordinates of pairs of variants
        for i in range(len(variants_pheno[gene]) -1):
            # if pair of variants not found, check the next variant
            if not_found:
                # grab the next variant
                for j in range(i+1, len(variants_pheno[gene])):
                    # check that variant do not overlap
                    var_coord1 = set(range(variants[variants_pheno[gene][i]][1], variants[variants_pheno[gene][i]][2]+1))
                    var_coord2 = set(range(variants[variants_pheno[gene][j]][1], variants[variants_pheno[gene][j]][2]+1))
                    if len(var_coord1.intersection(var_coord2)) == 0:
                        # keep both variants and exit loop
                        non_host_genes[gene] = [variants_pheno[gene][i], variants_pheno[gene][j]]
                        # variants are found
                        not_found = False
                        # exit inner loop
                        break
            else:
                # get the next gene
                break
    
    return non_host_genes
    
    
    
    
def variants_phenotype_host_genes(host_nested_variants_file, linkage_pickled, all_genes_phenotypes_file, gff_file):
    '''
    (file, file, file, file) -> dict
    Return a dict with host genes as key and a 2-item list as value containing
    2 non-overlapping phenotype-causing variants, one affecting a nested gene
    and one which doesn't      
    '''
    
    # get the coordinates of the variants
    # {var_ID: [chromo, start, end, source, type]}
    variants = variants_coord(gff_file)
    
    # get the variants contained within the host and affecting the host CDS
    # but without affecting the nested gene
    # { host_gene: [variant1, variant2]}
    hosts_variants_only = non_overlapping_variants(host_nested_variants_file, linkage_pickled, gff_file)
    
    # get the protein-coding genes with phenotype
    # create dict {gene: {var: [pheno1, pheno2]}}
    alleles = alleles_phenotypes(all_genes_phenotypes_file, linkage_pickled)
    
    # create a set of variants with phenotypes
    pheno_var = set()
    for gene in alleles:
        for var in alleles[gene]:
            pheno_var.add(var)
            
    # delete variants that have no phenotype
    for gene in hosts_variants_only:
        # create list of variant to delete
        to_remove = []
        # check if variant have a phenotype, if not then add to list
        for var in hosts_variants_only[gene]:
            if var not in pheno_var:
                to_remove.append(var)
        # delete variants without a phenotype
        for var in to_remove:
            hosts_variants_only[gene].remove(var)
    # remove genes without any variants
    to_remove = []
    for gene in hosts_variants_only:
        if len(hosts_variants_only[gene]) == 0:
            to_remove.append(gene)
    for gene in to_remove:
        del hosts_variants_only[gene]
        
    # get the variants with phenotype affecting the host and the nested gene
    # {host_gene: {nested_gene: [variant1, variant2]}}
    hosts_nested = host_nested_variant_relationship(host_nested_variants_file)
    
    # keep a single variant affecting the host and any nested gene
    # group all variant reagrdless of the nested gene they affect
    # {host_gene: [variant1, variant2]}
    hosts = {}
    for gene in hosts_nested:
        # initialize list
        hosts[gene] = []
        #loop over all nested genes and add variants to list
        for nested_gene in hosts_nested[gene]:
            hosts[gene].extend(hosts_nested[gene][nested_gene])
    
    # for each gene, find a pair of non-overlapping variants
    # create dict with host gene as key and 2-item list as value
    # containing a variant affecting the nested gene
    # and one not affeecting the nested gene
    host_genes = {}
    
    # loop over the genes in hosts
    for gene in hosts:
        # check that gene has variant not affecting the nested gene
        if gene in hosts_variants_only:
            # set up boolean
            not_found = True
            # compare coordinates of variants
            for i in range(len(hosts_variants_only[gene])):
                # if pair of variants not found, check the next variant
                if not_found:
                    for j in range(len(hosts[gene])):
                        # check that variants do not overlap
                        var_coord1 = set(range(variants[hosts_variants_only[gene][i]][1], variants[hosts_variants_only[gene][i]][2]+1))
                        var_coord2 = set(range(variants[hosts[gene][j]][1], variants[hosts[gene][j]][2]+1))
                        if len(var_coord1.intersection(var_coord2)) == 0:
                            # keep both variants and exit loop
                            host_genes[gene] = [hosts_variants_only[gene][i], hosts[gene][j]]
                            # variants are found
                            not_found = False
                            # exit inner loop
                            break
                else:
                    # get the next gene
                    break
                
    return host_genes
            
            
def jaccard_similarity(set1, set2):
    '''
    (set, set) -> float
    Take 2 sets and return the Jaccard similarity index, defined as the 
    number of common items in the 2 sets devided by the total number of
    uniue items in both sets
    '''
        
    jaccard = len(set1.intersection(set2)) / len(set1.union(set2))
    return jaccard


def phenotypes_of_variants(all_genes_phenotypes_file):
    '''
    (file) -> dict
    Return a dictionnary with variant ID as key and a set of associated 
    phenotype IDs
    '''
    
    # open file for reading
    infile = open(all_genes_phenotypes_file, 'r')
    # skip header
    infile.readline()
    
    # create dict {var: [pheno1, pheno2]}
    variants_phenotypes = {}
    
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
                if item.startswith('WBVar'):
                    var = item
                elif item.startswith('WBPheno'):
                    pheno = item
            # populate dict
            if var in variants_phenotypes:
                variants_phenotypes[var].add(pheno)
            else:
                variants_phenotypes[var] = {pheno}
                
    #close file after reading
    infile.close()
    
    return variants_phenotypes

                      
def phenotype_overlap(host_genes, non_host_genes, variants_phenotypes):
    '''
    (dict, dict, dict) -> (list, list)
    Take 3 dictionnaries: a dict for host-genes and a dict for non-host-genes,
    each with genes as keys and a 2-item list as value with variants affecting
    the CDS of the gene, and a dict of variants : phenotypes pairs.
    Return a tuple with lists of jaccard similarity indices respectively for 
    host-genes and non-host-genes
    '''
       
    # make a list of jaccard similarity indices for the host-genes
    host_genes_similarity = []
    # loop over host-genes
    for gene in host_genes:
        # make a set of phenotype for each variant
        var1 = variants_phenotypes[host_genes[gene][0]]
        var2 = variants_phenotypes[host_genes[gene][1]]
        # compute the jaccard similarity index
        jaccard = jaccard_similarity(var1, var2)
        # add index to list
        host_genes_similarity.append(jaccard)
    
    # make a list of jaccard similarity indices for the non-host-genes
    non_host_genes_similarity = []
    # loop over the non-host-genes
    for gene in non_host_genes:
        # make a set of phenotype for each variant
        var1 = variants_phenotypes[non_host_genes[gene][0]]
        var2 = variants_phenotypes[non_host_genes[gene][1]]
        # compute the jaccard similarity index
        jaccard = jaccard_similarity(var1, var2)
        # add index to list
        non_host_genes_similarity.append(jaccard)
    
    return host_genes_similarity, non_host_genes_similarity
        
        



                    
                    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    