# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 10:27:38 2015

@author: Richard
"""

import pickle



def all_elegans_genes(gff_file):
    '''
    (file) -> dict
    Return a dictionnary with chromosome as key and a dictionnary of gene
    attributes as value
    '''
      
    # open file to read
    infile = open(gff_file, 'r')

    # make a dictionnary with chromosome as key and a dictionnary as value
    # {'I': {gene:[chromosome, start, end, sense, biotype]}}
    linkage = {}
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[1] == 'WormBase' and 'ID=Gene' in line[8]:
                # get the gene name
                gene = line[8][line[8].index('ID=Gene:')+8: line[8].index(';')]
                # get the biotype
                biotype = line[8][line[8].index('biotype=')+8: line[8].index(';', line[8].index('biotype=')+8)]  
                # chromosome
                chromo = line[0]
                # start position
                start = int(line[3])
                # end position
                end = int(line[4])
                # sense
                sense = line[6]            
                # populate the dictionnary 
                if chromo not in linkage:
                    # create a dictionnary
                    linkage[chromo] = {}
                    linkage[chromo][gene] = [chromo, start, end, sense, biotype]
                else:
                    linkage[chromo][gene] = [chromo, start, end, sense, biotype]
    
    # close file after reading
    infile.close()
    
    return linkage


def find_nested_genes(linkage):
    '''
    (dict) -> dict
    Take a dictionnary of gene coordinates and return a dictionnary of host protein
    coding genes including nested genes
    '''
    
    # create a dict to store the host: nested gene pairs    
    host_genes = {}
    
    # loop over each chromosome
    for LG in linkage:
        print(LG)
        # make a list of genes located on LG
        Genes = [gene for gene in linkage[LG]]
        # compare each gene with all other genes and ask if genes are nested
        for i in range(len(Genes)):
            for j in range(i+1, len(Genes)):
                gene1, gene2 = Genes[i], Genes[j]
                # at least gene1 or gene2 needs to be a protein_coding gene
                if linkage[LG][gene1][-1] == 'protein_coding' or linkage[LG][gene2][-1] == 'protein_coding':
                    coordinates1 = set(range(linkage[LG][gene1][1], linkage[LG][gene1][2]+1))
                    coordinates2 = set(range(linkage[LG][gene2][1], linkage[LG][gene2][2]+1))
                    # compare gene coordinates
                    if coordinates1.issubset(coordinates2) and linkage[LG][gene2][-1] == 'protein_coding':
                        # gene1 is nested in gene2
                        if gene2 in host_genes:
                            host_genes[gene2].append(gene1)
                        else:
                            host_genes[gene2] = [gene1]
                    elif coordinates2.issubset(coordinates1) and linkage[LG][gene1][-1] == 'protein_coding':
                        # gene2 is nested in gene1
                        if gene1 in host_genes:
                            host_genes[gene1].append(gene2)
                        else:
                            host_genes[gene1] = [gene2]
                        
    return host_genes
 
    
def mutations(gff_file):
    '''
    (file) -> dict
    Return a dictionnary with chromosome as key and a dictionnary of variant
    features as value
    '''
    
    # open file for reading
    infile = open(gff_file, 'r')
    
    # create a dict to store the mutations {LG: {WBVAr: [start, end, variation_type]}}
    variants = {}
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
            # record only non single point mutations
            if start != end:
                if chromo not in variants:
                    variants[chromo] = {}
                    variants[chromo][var_name] = [start, end, source]
                else:
                    variants[chromo][var_name] = [start, end, source]
    # close file
    infile.close()
    
    return variants
                    

def gene_to_transcripts(gff_file):
    '''
    (file) -> dict
    Return a dictionnary with gene : transcripts pairs using the annotations
    of protein-coding genes in the gff file
    '''
        
    # open file for reading
    infile = open(gff_file, 'r')
    # create a dictionnary 
    genes = {}
    # loop over file, set gene as key and a list of transcript as value
    for line in infile:
        # use only WormBase valid genes
        if 'WormBase' in line:
            line = line.rstrip().split()
            # check that source is WormBase
            if line[1] == 'WormBase':
                # grab the gene name
                if line[2] == 'gene' and 'ID=Gene' in line[8] and 'protein_coding' in line[8]:
                    gene = line[8][line[8].index('ID=Gene:')+8: line[8].index(';')]
                    genes[gene] = []
                # grab the transcripts
                elif line[2] == 'mRNA':
                    transcript = line[8][line[8].index('ID=Transcript:')+14: line[8].index(';')]
                    gene = line[8][line[8].index('Gene:')+5 : line[8].index(';', line[8].index('Gene:'))]
                    if gene in genes:
                        genes[gene].append(transcript)
    #close file after reading                
    infile.close()
    return genes


def transcripts_to_gene(gff_file):
    '''
    (file) -> dict
    Return a dictionnary of transcript : gene pairs using the annotations of protein-coding genes in the gff file
    '''
    # create a dictionnary of gene : list of transcripts pairs   
    genes = gene_to_transcripts(gff_file)
    
    # reverse the genes dict
    transcripts = {}
    for gene, ts in genes.items():
        for i in range(len(ts)):
            transcripts[ts[i]] = gene
            
    return transcripts
    

def coding_sequences(gff_file):
    '''
    (file) -> dict
    Return a dictionnary with chromosome as key and a dictionnary of gene name
    and CDS positions list pairs as value    
    '''
    
    # create a dictionnary of transcripts: gene pairs
    transcripts = transcripts_to_gene(gff_file)
    
    # open file for reading
    infile = open(gff_file, 'r')
        
    # make a dictionnary with chromosome as key and a dictionnary as value
    # {'I': {gene:[(CDS_start, CDS_end), (CDS_start, CDS_end)]}}
    chromosomes = {}
    
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[1] == 'WormBase':
                if line[2] == 'CDS':
                    # get the transcript name of that CDS
                    # chop line[8] up to first transcript
                    line[8] = line[8][line[8].index('Transcript'):]
                    # split on ';' and remove any str not containing Transcript
                    line[8] = line[8].split(';')
                    to_remove = []
                    for item in line[8]:
                        if 'Transcript' not in item:
                            to_remove.append(item)
                    for item in to_remove:
                        line[8].remove(item)
                    # get remaining string
                    ts_name = line[8][0]
                    # remove 'Transcript:' from names
                    ts_name = ts_name.replace('Transcript:', '')
                    # check is names contain multiple transcript names
                    if ',' in ts_name:
                        # split names on ',':
                        ts_name = ts_name.split(',')
                        # we only need to record all the CDS, doesn't matter
                        # if multiple transcripts of the same gene share the same CDS
                        ts_name = ts_name[0]
                    # get the start and end positions
                    start, end = int(line[3]), int(line[4])
                    # get the chromosome
                    chromo = line[0]
                    # get the parent gene of transcript
                    gene = transcripts[ts_name]
                    # check if chromo in chromosomes
                    if chromo in chromosomes:
                        # check if parent gene in chromosomes[chromo]
                        if gene in chromosomes[chromo]:
                            # add the CDS coordinates
                            chromosomes[chromo][gene].append([start, end])
                        else:
                            # initialize list
                            chromosomes[chromo][gene] = [[start, end]]
                    else:
                        # initialize dict
                        chromosomes[chromo] = {}
                        chromosomes[chromo][gene] = [[start, end]]
    # close file after reading
    infile.close()
    
    return chromosomes



def CDS_host_genes(coding, host_genes):
    '''
    (dict, dit) -> dict
    Take a dictionnary with CDS coordinates for all genes and a dictionnary
    of host genes including nested genes and return a dictionnary with the CDS
    coordinates of the host genes only
    '''
    
    # create a restricted dict of coding sequences for host genes
    host_genes_coding = {}
    # loop over chromosomes in coding
    for chromo in coding:
        # loop over all genes on chromo,
        # if gene is host then add coding seq to host_genes_coding
        for gene in coding[chromo]:
            if gene in host_genes:
                # initiate dict if chromo not in host_genes_coding
                if chromo in host_genes_coding:
                    # make a list with the coding[LG][gene]
                    host_genes_coding[chromo][gene] = list(coding[chromo][gene])
                else:
                    # create a chromo : dict pair
                    host_genes_coding[chromo] = {}
                    # make a list with the coding[LG][gene]
                    host_genes_coding[chromo][gene] = list(coding[chromo][gene])
    
    return host_genes_coding
    

def variants_host_coding(variants, host_genes_coding):
    '''
    (dict, dict) -> dict
    Take a dictionnary with the coordinates of non-single point mutation
    variants, and a dictionnary including the coordinates of the host genes
    CDS regions and return a dictionnary with variants affecting the coding
    regions of host genes   
    '''
    # make a restricted dict with variants affecting the CDS of host genes
    # {chromo: {gene: [(variant1, start, end, source), (variant2, start, end, source)]}}
    host_CDS_variants = {}
    # loop over chromosome in variants:
    for chromo in variants:
        print(chromo)
        # loop over all variants on chromo
        for variant_name in variants[chromo]:
            # get variant coordinates
            variant_coord = set(range(variants[chromo][variant_name][0], variants[chromo][variant_name][1]+1))
            # compare variant positions to CDS positions for each gene on chromo
            # make sure that host protein_coding genes are on chromo
            if chromo in host_genes_coding:
                for gene in host_genes_coding[chromo]:
                    # get the coordinates of each CDS
                    for CDS in host_genes_coding[chromo][gene]:
                        CDS_coord = set(range(CDS[0], CDS[1]+1))
                        # check if gene's CDS overlap with variant coordinates
                        if len(variant_coord.intersection(CDS_coord)) != 0:
                            # if overlap, create a list [variant_name, start, end, source]
                            variant_info = [variant_name]
                            for i in range(len(variants[chromo][variant_name])):
                                variant_info.append(variants[chromo][variant_name][i])
                            # check if chromo in dict
                            if chromo in host_CDS_variants:
                                # check if gene is key of host_CDS_variants[chromo]
                                if gene in host_CDS_variants[chromo]:
                                    # add [variant_name, start, end, source] to list
                                    host_CDS_variants[chromo][gene].append(variant_info)
                                else:
                                    # initialize list
                                    host_CDS_variants[chromo][gene] = [variant_info]
                            else:
                                # initialize dict
                                host_CDS_variants[chromo] = {}
                                # initialize list
                                host_CDS_variants[chromo][gene] = [variant_info]
                            # no need to evaluate all CDS, it needs only 1 CDS to be affected by the variant
                            break
    
    return host_CDS_variants


def host_nested_alterations(host_CDS_variants, host_genes, linkage):
    '''
    (dict, dict, dict) -> dict
    Take a dictionnary of variants affecting the CDS of host genes, a 
    dictionnary of host: nested gene pairs and a dictionnary with all the genes
    coordinates and return a dictionnary of host genes with dictionnaries of
    nested genes and list of variant pairs as value. Variants are alterations
    of the CDS of the host gene and alterations of the nested gene
    '''
        
    # make a dictionnary of variants affecting the host gene's CDS and the nested gene
    # { host_gene: {nested_gene: [[variant1, start, end, source], [variant2, start, end, source]]}}
    host_nested_variants = {}
    
    # loop over chromosomes in host_CDS_variants
    for chromo in host_CDS_variants:
        # loop over all gene on chromo
        for gene in host_CDS_variants[chromo]:
            # check all variants affecting the CDS of that gene
            for variant in host_CDS_variants[chromo][gene]:
                # get the variant coordinates
                variant_coord = set(range(variant[1], variant[2]+1))
                # compare the variant coordinates with the nested gene coordinates
                for nested_gene in host_genes[gene]:
                    # get the nested genes coordinates
                    nested_coord = set(range(linkage[chromo][nested_gene][1], linkage[chromo][nested_gene][2]+1))
                    # check if variant and nested gene coordinates overlap
                    if len(variant_coord.intersection(nested_coord)) != 0:
                        # check if gene in host_nested_variants
                        if gene in host_nested_variants:
                            # check if nested_genes is key
                            if nested_gene in host_nested_variants[gene]:
                                host_nested_variants[gene][nested_gene].append(variant)
                            else:
                                host_nested_variants[gene][nested_gene] = [variant]
                        else:
                            # initialize dict
                            host_nested_variants[gene] = {}
                            host_nested_variants[gene][nested_gene] = [variant]
    
    return host_nested_variants

                    
def host_nested_variants_to_file(gff_file, outputfile):
    '''
    (file) -> file
    Find the variants altering both the host genes CDS and the sequence of the
    nested gene and save information about the host, nested genes and variants
    to an outputfile
    '''
       
    # get the positions of all genes
    linkage = all_elegans_genes(gff_file)
    print('collected coordinates of all annotated genes')
    
    # find the nested and host genes genes
    host_genes = find_nested_genes(linkage)
    print('found the host-nested pairs')
    
    # get the coordinates of the coding sequence of all protein-coding genes
    coding = coding_sequences(gff_file)
    print('got the genes\' coordinates')
    
    # get the CDS coordinates of the host genes
    host_genes_coding = CDS_host_genes(coding, host_genes)
    print('got the host genes\' CDS coordinates')
    
    # get variant coordinates
    variants = mutations(gff_file)
    print('collected all non-SNP variants')
    
    # make a dictionnary with variants affecting the host genes' CDS
    host_CDS_variants = variants_host_coding(variants, host_genes_coding)
    print('found variants affecting the host genes\' CDS')
    
    # make a dict with variants affecting the host gene CDS and the nested gene
    host_nested_variants = host_nested_alterations(host_CDS_variants, host_genes, linkage)
    print('collected variants affecting host and nested genes')
    
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
                        nested_start = str(linkage[LG][nested_gene[1]])
                        nested_end = str(linkage[LG][nested_gene[2]])
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
                variation += (variant + '_')
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
    

def nested_genes_per_host(host_genes_pickled):
    '''
    (file) -> dict
    Return a dict with the number of nested genes per host as key and the 
    host genes count as value    
    '''
    
    # unpickle the host gene dict
    hostfile = open(host_genes_pickled, 'rb')
    host_genes = pickle.load(hostfile)
    hostfile.close()
    
    # create dict of # nested genes : count pairs
    nested = {}
    for gene in host_genes:
        num_nested = len(host_genes[gene])
        if num_nested in nested:
            nested[num_nested] += 1
        else:
            nested[num_nested] = 1
    
    return nested


def unpickle_dict(pickled_dict):
    '''
    (file) -> dict
    Return a dict by unpickling the dict stored in pickled file
    '''
    
    # unpickle dicts
    infile = open(pickled_dict, 'rb')
    dictionnary = pickle.load(infile)
    infile.close()
    
    return dictionnary



def nested_host_pairs(host_genes_pickled):
    '''
    (file) -> dict
    Reverse the host_genes dictionnary and return a dictionnary of 
    nested gene : host gene pairs
    '''
    
    # unpickle dict with host_genes
    host_genes = unpickle_dict(host_genes_pickled)
    
    # create a dict of nested_gene : [host_genes] pairs
    nested_genes = {}
    for host in host_genes:
        for gene in host_genes[host]:
            if gene in nested_genes:
                nested_genes[gene].append(host)
            else:
                nested_genes[gene] = [host]
    
    return nested_genes


def gene_coordinates(linkage_pickled):
    '''
    (file) -> dict
    Return a dictionnary with gene as key and coordinates as value
    '''
    
    # unpickle dict
    linkage = unpickle_dict(linkage_pickled)
        
    # create dict of gene coordinates
    # {gene:[chromosome, start, end, sense, biotype]}
    all_genes = {}
    for chromo in linkage:
        for gene in linkage[chromo]:
            all_genes[gene] = list(linkage[chromo][gene])

    return all_genes


def mutually_exclusive_gene_sets(host_genes_pickled, linkage_pickled):
    '''
    (file, file) -> (set, set, set)
    Return a tuple containing 3 sets of mutually exclusive protein coding genes
    that are hosts, nested and neither of these two categories
    '''
    
    # unpickle dict
    host_genes = unpickle_dict(host_genes_pickled)
    # create a dict of nested gene : host gene pairs
    nested_genes = nested_host_pairs(host_genes_pickled)
    # get the genes' coordinates
    all_genes = gene_coordinates(linkage_pickled)    
    
    # create sets of mutually exclusive genes
    nested_only = set()
    host_only = set()
    other_genes = set()
    
    # loop over host genes    
    for gene in host_genes:
        # check that gene is protein coding
        if all_genes[gene][-1] == 'protein_coding':
            # check that gene is not nested
            if gene not in nested_genes:
                # add gene to set of hosts
                host_only.add(gene)
                
    # loop over nested genes
    for gene in nested_genes:
        # verify that gene is protein coding
        if all_genes[gene][-1] == 'protein_coding':
            # verify that gene is not host
            if gene not in host_genes:
                # add gene to set of nested genes
                nested_only.add(gene)
                
    # loop over all genes
    for gene in all_genes:
        # verify that gene is protein coding
        if all_genes[gene][-1] == 'protein_coding':
            # verify that gene is not host and not nested
            if gene not in nested_genes and gene not in host_genes:
                # add gene to set of other genes
                other_genes.add(gene)
                
    return host_only, nested_only, other_genes
                
        
def host_nested_others_gene_length(host_genes_pickled, linkage_pickled):
    '''
    (file, file) -> (list, list, list)
    Return a tuple containing lists of gene length for hosts, nested genes and
    genes that are neither host or nested genes
    '''
    
    # create a sets of protein-coding genes that are mutually exclusive
    host_only, nested_only, other_genes = mutually_exclusive_gene_sets(host_genes_pickled, linkage_pickled)
    # get the genes' coordinates
    all_genes = gene_coordinates(linkage_pickled) 
    # create a list to store the length of host genes
    host_length = [(all_genes[gene][2] - all_genes[gene][1] + 1) for gene in host_only]
    # create a list to store the length of nested genes
    nested_length = [(all_genes[gene][2] - all_genes[gene][1] + 1) for gene in nested_only]
    # create a list to store the length of non-host and non-nested genes
    others_length = [(all_genes[gene][2] - all_genes[gene][1] + 1) for gene in other_genes]
    
    return host_length, nested_length, others_length
    

def host_nested_others_exon_length(host_genes_pickled, linkage_pickled, gff_file):
    '''
    (file, file) -> (list, list, list)
    Return a tuple containing a list of exon length for hosts, nested genes and
    genes that are neither host or nested
    '''
    
    # create a sets of protein-coding genes that are mutually exclusive
    host_only, nested_only, other_genes = mutually_exclusive_gene_sets(host_genes_pickled, linkage_pickled)    
    
    # get the CDS coordinates
    # {chromo: {gene:[(CDS_start, CDS_end), (CDS_start, CDS_end)]}}
    coding = coding_sequences(gff_file)
    
    # create lists to store the CDS length
    host_length = []
    nested_length = []
    others_length = []
        
    # loop over chromosomes
    for chromo in coding:
        # loop over genes
        for gene in coding[chromo]:
            # compute CDS length
            for CDS in coding[chromo][gene]:
                CDS_length = CDS[1] - CDS[0] + 1
                # check if gene is host or not, add CDS length to list
                if gene in host_only:
                    host_length.append(CDS_length)
                elif gene in nested_only:
                    nested_length.append(CDS_length)
                elif gene in other_genes:
                    others_length.append(CDS_length)
                                  
    return host_length, nested_length, others_length            
                
    



#####################
    
#### ADD an option to coding_sequences function to record intron sequences instead
#### of CDS. modify other functions
    
def intron_sequences(gff_file):
    '''
    (file) -> dict
    Return a dictionnary with chromosome as key and a dictionnary of gene name
    and intron positions list pairs as value    
    '''
    
    # create a dictionnary of transcripts: gene pairs
    transcripts = transcripts_to_gene(gff_file)
    
    # open file for reading
    infile = open(gff_file, 'r')
        
    # make a dictionnary with chromosome as key and a dictionnary as value
    # {'I': {gene:[(intron_start, intron_end), (intron_start, intron_end)]}}
    chromosomes = {}
    
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[1] == 'WormBase':
                if line[2] == 'intron':
                    # get the transcript name of that CDS
                    # chop line[8] up to first transcript
                    line[8] = line[8][line[8].index('Transcript'):]
                    # split on ';' and remove any str not containing Transcript
                    line[8] = line[8].split(';')
                    to_remove = []
                    for item in line[8]:
                        if 'Transcript' not in item:
                            to_remove.append(item)
                    for item in to_remove:
                        line[8].remove(item)
                    # get remaining string
                    ts_name = line[8][0]
                    # remove 'Transcript:' from names
                    ts_name = ts_name.replace('Transcript:', '')
                    # check is names contain multiple transcript names
                    if ',' in ts_name:
                        # split names on ',':
                        ts_name = ts_name.split(',')
                        # we only need to record all the CDS, doesn't matter
                        # if multiple transcripts of the same gene share the same CDS
                        ts_name = ts_name[0]
                    # get the start and end positions
                    start, end = int(line[3]), int(line[4])
                    # get the chromosome
                    chromo = line[0]
                    # get the parent gene of transcript
                    gene = transcripts[ts_name]
                    # check if chromo in chromosomes
                    if chromo in chromosomes:
                        # check if parent gene in chromosomes[chromo]
                        if gene in chromosomes[chromo]:
                            # add the CDS coordinates
                            chromosomes[chromo][gene].append([start, end])
                        else:
                            # initialize list
                            chromosomes[chromo][gene] = [[start, end]]
                    else:
                        # initialize dict
                        chromosomes[chromo] = {}
                        chromosomes[chromo][gene] = [[start, end]]
    # close file after reading
    infile.close()
    
    return chromosomes  
    

def host_nested_others_intron_counts(host_genes_pickled, linkage_pickled, gff_file):
    '''
    (file, file) -> (list, list, list)
    Return a tuple containing lists of intron counts for host genes, nested genes
    and genes that are neither host or nested
    '''
    
    # create a sets of protein-coding genes that are mutually exclusive
    host_only, nested_only, other_genes = mutually_exclusive_gene_sets(host_genes_pickled, linkage_pickled) 
    
    # get the intron coordinates
    # {chromo: {gene:[(intron_start, intron_end), (intron_start, intron_end)]}}
    introns = intron_sequences(gff_file)
    
    # create a dict {gene: [(intron_start, intron_end), (intron_start, intron_end)]}
    genes_with_introns = {}
    for chromo in introns:
        for gene in introns[chromo]:
            genes_with_introns[gene] = list(introns[chromo][gene])
    
    # create lists to store the number of introns
    host_introns = []
    nested_introns = []
    others_introns = []
    
    # loop over hosts
    for gene in host_only:
        # if gene has intron, add the number of introns to list
        if gene in genes_with_introns:
            host_introns.append(len(genes_with_introns[gene]))
        else:
            # gene is intronless, add 0
            host_introns.append(0)
    
    # loop over nested genes
    for gene in nested_only:
        # if gene has intron, add the number of introns to list
        if gene in genes_with_introns:
            nested_introns.append(len(genes_with_introns[gene]))
        else:
            # gene is intronless, add 0
            nested_introns.append(0)
            
    # loop over other genes
    for gene in other_genes:
        # if gene has intron, add the number of introns to list
        if gene in genes_with_introns:
            others_introns.append(len(genes_with_introns[gene]))
        else:
            # gene is intronless, add 0
            others_introns.append(0)
    
    return host_introns, nested_introns, others_introns
    
    
     

def host_nested_others_intron_length(host_genes_pickled, linkage_pickled, gff_file):
    '''
    (file, file) -> (list, list, list)
    Return a tuple containing lists of intron lengths for host-genes, nested 
    genes and genes that are not host and nested.
    Precondition, do not consider intronless genes
    '''
    
    # create a sets of protein-coding genes that are mutually exclusive
    host_only, nested_only, other_genes = mutually_exclusive_gene_sets(host_genes_pickled, linkage_pickled)     
    
    # get the introns coordinates
    # {chromo: {gene:[(intron_start, intron_end), (intron_start, intron_end)]}}
    introns = intron_sequences(gff_file)
    
    # create lists to store the intron length
    host_introns = []
    nested_introns = []
    others_introns = []    
    
    # loop over chromosomes
    for chromo in introns:
        # loop over genes
        for gene in introns[chromo]:
            # compute intron length
            for intron in introns[chromo][gene]:
                intron_length = intron[1] - intron[0] + 1
                # check if gene is host or not, add intron length to list
                if gene in host_only:
                    host_introns.append(intron_length)
                elif gene in nested_only:
                    nested_introns.append(intron_length)
                elif gene in other_genes:
                    others_introns.append(intron_length)
                                  
    return host_introns, nested_introns, others_introns
    
  

def ratio_exon_intron(host_genes_pickled, linkage_pickled, gff_file):
    '''
    (file, file, file) -> (list, list, list)
    Return a tuple containing lists of ratio of total intron length over total 
    exon length for host genes, nested genes and genes that are not host or nested
    '''
    
    # create a sets of protein-coding genes that are mutually exclusive
    host_only, nested_only, other_genes = mutually_exclusive_gene_sets(host_genes_pickled, linkage_pickled)     
    
    # create lists to store the ratio of inton length over exon length
    host_ratio = []
    nested_ratio = []
    others_ratio = []    

    # get the introns coordinates
    # {chromo: {gene:[(intron_start, intron_end), (intron_start, intron_end)]}}
    introns = intron_sequences(gff_file)
     
    # create a dict {gene: [(intron_start, intron_end), (intron_start, intron_end)]}
    genes_with_introns = {}
    for chromo in introns:
        for gene in introns[chromo]:
            genes_with_introns[gene] = list(introns[chromo][gene])
            
    # get the CDS coordinates
    # {chromo: {gene:[(CDS_start, CDS_end), (CDS_start, CDS_end)]}}
    coding = coding_sequences(gff_file)
    
    # loop over chromosomes in coding
    for chromo in coding:
        # loop over genes on chromosome
        for gene in coding[chromo]:
            # compute total CDS length
            total_CDS_length = 0
            for CDS in coding[chromo][gene]:
                CDS_length = CDS[1] - CDS[0] + 1
                total_CDS_length += CDS_length
            # compute total intron length
            total_intron_length = 0
            # check if gene is intronless or not
            if gene in genes_with_introns:
                for intron in genes_with_introns[gene]:
                    intron_length = intron[1] - intron[0] + 1
                    total_intron_length += intron_length
            # compute ratio
            ratio = total_intron_length / total_CDS_length
            # check if gene is host nested or other, add ratio to list
            if gene in host_only:
                host_ratio.append(ratio)
            elif gene in nested_only:
                nested_ratio.append(ratio)
            elif gene in other_genes:
                others_ratio.append(ratio)
                
    return host_ratio, nested_ratio, others_ratio


    
def nested_genes_location(host_genes_pickled, linkage_pickled, gff_file):
    '''
    (file) -> tuple
    Return a tuple with counts of nested genes located within introns, 
    within coding-exons, or other location within the host genes
    (ex: overlapping exons and introns, UTRs, non-coding exons)
    '''
    
    # create dict of gene coordinates
    # {gene:[chromosome, start, end, sense, biotype]}
    all_genes = gene_coordinates(linkage_pickled)
        
    # create a dict of nested_gene : host_gene pairs
    nested_genes = nested_host_pairs(host_genes_pickled)
        
    # get the coordinates of introns for all protein-coding genes
    # {chromo: {gene: [[start, end], [start, end]]}}
    introns = intron_sequences(gff_file)
    
    # get the coordinates of CDS for all protein-coding genes
    # {chromo: {gene: [[start, end], [start, end]]}}    
    coding = coding_sequences(gff_file)
    
    # create set of nested genes located in introns
    intronic = set()
    
    # create set of nested genes located in CDS
    exonic = set()
    
    # create set of nested genes not intronic and not exonic
    others = set()
    
    # loop over nested genes
    for gene in nested_genes:
        # get the host gene
        host = nested_genes[gene]
        # get chromosome of host and nested
        chromo = all_genes[gene][0]
        # get the nested genes' coordinates
        nested_start, nested_end = all_genes[gene][1], all_genes[gene][2]
        nested_coord = set(range(nested_start, nested_end + 1))
        # compare the nested gene coordinates with the intron coordinates
        if host in introns[chromo]:
            for intron in introns[chromo][host]:
                # get the introns coordinates
                intron_start, intron_end = intron[0], intron[1]
                intron_coord = set(range(intron_start, intron_end + 1))
                # check if nested gene is entirely contained within intron
                if nested_coord.issubset(intron_coord):
                    intronic.add(gene)
                    # exit inner loop if intronic nested gene
                    break
        # compare the nested gene coordinates with the CDS coordinates
        if host in coding[chromo]:
            for CDS in coding[chromo][host]:
                # get the CDS coordinates
                CDS_start, CDS_end = CDS[0], CDS[1]
                CDS_coord = set(range(CDS_start, CDS_end + 1))
                # check if nested gene is entirely contained within CDS
                # check if nested gene not already in intronic
                # 19 genes can be exonic or intronic because of alternative splicing
                if nested_coord.issubset(CDS_coord) and gene not in intronic:
                    exonic.add(gene)
                    # exit inner loop if exonic nested gene
                    break
        # classify as other if nested gene is not exonic and not intronic
        if gene not in intronic and gene not in exonic:
            others.add(gene)
            
    return intronic, exonic, others
                
    
    
    
    
    