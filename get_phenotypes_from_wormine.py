# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 14:18:42 2015

@author: Richard
"""

# this code is used to get the phenotypes for alleles of the host genes
# by using Wormine API intermine

from intermine.webservice import Service

# fetch data
service = Service("http://www.wormbase.org/tools/wormmine/service")

newfile = open('Celegans_alleles_phenotypes.txt', 'w')

header = '\t'.join(['WormBase_Gene_ID', 'Sequence_Name', 'Gene_Name', 'Alleles_WormBase_ID', 
                    'Alleles_Public_Name', 'Phenotypes_Observed_Identifier', 'Phenotypes_Observed_Name']) 
newfile.write(header + '\n')

# create a set of genes that have a phenotype
pheno = open('genes_with_phenotypes.txt', 'r')
# create a set of host genes
pheno_genes = set()
for line in pheno:
    line = line.rstrip()
    if line != '':
        pheno_genes.add(line)


# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("Gene.primaryIdentifier", "ASC")

for gene in pheno_genes:
    
    # Get a new query on the class (table) you will be querying:
    query = service.new_query("Gene")
    
    
    # The view specifies the output columns
    query.add_view(
        "primaryIdentifier", "secondaryIdentifier", "symbol",
        "alleles.primaryIdentifier", "alleles.symbol",
        "alleles.phenotypesObserved.identifier", "alleles.phenotypesObserved.name"
    )
        
    # You can edit the constraint values below
    query.add_constraint("symbol", "=", gene, "A")

    # Uncomment and edit the code below to specify your own custom logic:
    # query.set_logic("A")
    
    for row in query.rows():
        newfile.write(row["primaryIdentifier"])
        newfile.write('\t')
        newfile.write(row["secondaryIdentifier"])
        newfile.write('\t')
        newfile.write(row["symbol"])
        newfile.write('\t')
        newfile.write(row["alleles.primaryIdentifier"])
        newfile.write('\t')
        newfile.write(row["alleles.symbol"])
        newfile.write('\t')
        newfile.write(row["alleles.phenotypesObserved.identifier"])
        newfile.write('\t')
        newfile.write(row["alleles.phenotypesObserved.name"])
        newfile.write('\n')
    
newfile.close()
        
